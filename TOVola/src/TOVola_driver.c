#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "gsl/gsl_odeiv2.h"
#include "gsl/gsl_errno.h"

#include "GRHayLib.h"
#include "TOVola_interp.h"
#include "TOVola_defines.h"
#include "TOVola_solve.h"

/********************************
// TOVola is a TOV solver inspired by NRPy designed specifically for use in the
Einstein Toolkit.
// It integrates the TOV equations for a solution to static, spherically
symmetric stars (Neutron stars in specific), and interpollates the solution the
the ET grid.
//
// MOST IMPORTANTLY, it handles multiple types of EOS in the solution: Simple
Polytrope, Piecewise Polytrope, and Tabulated EOS. An upgrade from the original
ET TOVSolver thorn.
//
// Simple Polytrope:
//      -- The Simple Polytrope is a special version of the Piecewise Polytrope
with one EOS region relating Pressure to baryon density (P=K*rho_baryon^Gamma),
where K and Gamma are parameters to be input
//      -- Only one K and Gamma are necessary to run this implementation, and a
Gamma_Thermal defaulted to 2.0
//      -- Hamiltonian Constraint Violation Validated with Baikal
// Piecewise Polytrope:
//      -- The more general verion of the simple polytrope where you can choose
the number of regions for your polytrope, each region denoted by
(P=K_i*rho_baryon^Gamma_i)
//      -- You can get your necessary parameters by following the method
outlined in Read et. al. (2008)
//      -- Since we follow Read's method, you only need to input one K value
(The GRHayL library will handle the rest), and all the Gammas and rho_baryon
boundary points for the number of specified regions.
//      -- Hamiltonian Constraint Violation Validated with Baikal
// Tabulated EOS:
//      -- The tabulated solver reads in a EOS table from your computer.
//      -- Once located, GRHayL slices the table for beta_equilibirum and then
uses an interpollator to find the necessary densities through the calculations.
//      -- Hamiltonian Constraint Violation Validated with Baikal
//
// Users who wish to interface with this program need only create a parfile for
a simulation. Output variables will be adjusted for use in ADMbase and Hydrobase
as initial data.
// IllinoisGRMHD is used to set Tmunu from the data and Baikal tests the
constraint Violation.
//
// Gallery examples are currently in the works, using GRHydro for evolution.
However, we have also used TOVola with GRHayL evolutions for a paper currently
in prep.
//
// David Boyer (10/08/24)
********************************/

#define ODE_SOLVER_DIM 4
#define TOVOLA_PRESSURE 0
#define TOVOLA_NU 1
#define TOVOLA_MASS 2
#define TOVOLA_R_ISO 3
#define NEGATIVE_R_INTERP_BUFFER 11

// Check to make sure the parameters given are valid.
void TOVola_Parameter_Checker(CCTK_ARGUMENTS) {

  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS_TOVola_Parameter_Checker;

  CCTK_INFO("TOVola Validating Interpolation stencil and EOS_type...");

  if (TOVola_Interpolation_Stencil > TOVola_Max_Interpolation_Stencil) {
    CCTK_PARAMWARN("TOVola_Interpolation_Stencil must not exceed the "
                   "Max_Interpolation_Stencil");
  }
  if (CCTK_EQUALS("Simple", TOVola_EOS_type) && ghl_eos->neos != 1) {
    CCTK_PARAMWARN("Error: Too many regions for the simple polytrope. Check "
                   "your value for neos, or use a piecewise polytrope.");
  }
  CCTK_INFO("Validation Complete!");
}

// Drive the TOV integration using GSL and then interpolate the data to the ET
// grid.
void TOVola_Solve_and_Interp(CCTK_ARGUMENTS) {

  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS_TOVola_Solve_and_Interp;

  TOVola_ID_persist_struct
      TOVola_ID_persist_tmp; // allocates memory for the pointer below.
  TOVola_ID_persist_struct *restrict TOVola_ID_persist = &TOVola_ID_persist_tmp;
  CCTK_REAL current_position = 0;

  /* Set up ODE system and driver */
  TOVola_data_struct TOVdata_tmp; // allocates memory for the pointer below.
  TOVola_data_struct *restrict TOVdata = &TOVdata_tmp;
  gsl_odeiv2_system system;
  gsl_odeiv2_driver *driver;

  // Setting EOS
  if (CCTK_EQUALS("Simple", TOVola_EOS_type)) {
    TOVdata->eos_type = 0;
  } else if (CCTK_EQUALS("Piecewise", TOVola_EOS_type)) {
    TOVdata->eos_type = 1;
  } else if (CCTK_EQUALS("Tabulated", TOVola_EOS_type)) {
    TOVdata->eos_type = 2;
    ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(TOVola_Tin, ghl_eos);
  }

  // Initialize other TOVdata member variables
  TOVdata->numpoints_actually_saved = 0;
  TOVdata->error_limit = TOVola_error_limit;
  TOVdata->initial_ode_step_size = TOVola_initial_ode_step_size;
  TOVdata->absolute_max_step = TOVola_absolute_max_step;
  TOVdata->absolute_min_step = TOVola_absolute_min_step;
  TOVdata->central_baryon_density = TOVola_central_baryon_density;
  TOVdata->ghl_eos = ghl_eos;
  if (setup_ode_system(TOVola_ODE_method, &system, &driver, TOVdata) != 0) {
    CCTK_ERROR("Failed to set up ODE system.");
  }

  /* Initialize ODE variables */
  CCTK_REAL TOVola_eq[ODE_SOLVER_DIM];
  CCTK_REAL c[2];
  TOVola_get_initial_condition(TOVola_eq, TOVdata);
  TOVola_assign_constants(c, TOVdata);

  /* Initial memory allocation */
  TOVdata->numels_alloced_TOV_arr = 1024;
  if (initialize_tovola_data(TOVdata) != 0) {
    gsl_odeiv2_driver_free(driver);
    CCTK_ERROR("Failed to initialize TOVola_data_struct.");
  }

  /* Integration loop */
  TOVdata->r_lengthscale =
      TOVola_initial_ode_step_size; // initialize dr to a crazy small value in
                                    // double precision.
  for (CCTK_INT i = 0; i < TOVola_size; i++) {
    CCTK_REAL dr = 0.01 * TOVdata->r_lengthscale;
    if (TOVdata->rho_baryon < 0.05 * TOVola_central_baryon_density) {
      // To get a super-accurate mass, reduce the dr sampling near the surface
      // of the star.
      dr = 1e-6 * TOVdata->r_lengthscale;
    }
    /* Exception handling */
    TOVola_exception_handler(current_position, TOVola_eq);

    /* Apply ODE step */
    CCTK_INT status = gsl_odeiv2_driver_apply(driver, &current_position,
                                              current_position + dr, TOVola_eq);
    if (status != GSL_SUCCESS) {
      CCTK_VINFO("GSL ODE solver failed with status %d.", status);
      gsl_odeiv2_driver_free(driver);
      CCTK_ERROR("Shutting down due to error...");
    };

    /* Post-step exception handling */
    TOVola_exception_handler(current_position, TOVola_eq);

    /* Evaluate densities */
    TOVola_evaluate_rho_and_eps(current_position, TOVola_eq, TOVdata);
    TOVola_assign_constants(c, TOVdata);

    /* Check if reallocation is needed */
    if (TOVdata->numpoints_actually_saved >= TOVdata->numels_alloced_TOV_arr) {
      // Update arr_size instead of modifying the macro
      const CCTK_INT new_arr_size = 1.5 * TOVdata->numels_alloced_TOV_arr;
      TOVdata->numels_alloced_TOV_arr = new_arr_size;
      TOVdata->rSchw_arr =
          realloc(TOVdata->rSchw_arr, sizeof(CCTK_REAL) * new_arr_size);
      TOVdata->rho_energy_arr =
          realloc(TOVdata->rho_energy_arr, sizeof(CCTK_REAL) * new_arr_size);
      TOVdata->rho_baryon_arr =
          realloc(TOVdata->rho_baryon_arr, sizeof(CCTK_REAL) * new_arr_size);
      TOVdata->P_arr =
          realloc(TOVdata->P_arr, sizeof(CCTK_REAL) * new_arr_size);
      TOVdata->M_arr =
          realloc(TOVdata->M_arr, sizeof(CCTK_REAL) * new_arr_size);
      TOVdata->nu_arr =
          realloc(TOVdata->nu_arr, sizeof(CCTK_REAL) * new_arr_size);
      TOVdata->Iso_r_arr =
          realloc(TOVdata->Iso_r_arr, sizeof(CCTK_REAL) * new_arr_size);

      if (!TOVdata->rSchw_arr || !TOVdata->rho_energy_arr ||
          !TOVdata->rho_baryon_arr || !TOVdata->P_arr || !TOVdata->M_arr ||
          !TOVdata->nu_arr || !TOVdata->Iso_r_arr) {
        CCTK_ERROR("Memory reallocation failed during integration.\n");
        gsl_odeiv2_driver_free(driver);
      }
    }

    /* Store data */
    TOVdata->rSchw_arr[TOVdata->numpoints_actually_saved] = current_position;
    TOVdata->rho_energy_arr[TOVdata->numpoints_actually_saved] = c[0];
    TOVdata->rho_baryon_arr[TOVdata->numpoints_actually_saved] = c[1];
    TOVdata->P_arr[TOVdata->numpoints_actually_saved] =
        TOVola_eq[TOVOLA_PRESSURE];
    TOVdata->M_arr[TOVdata->numpoints_actually_saved] = TOVola_eq[TOVOLA_MASS];
    TOVdata->nu_arr[TOVdata->numpoints_actually_saved] = TOVola_eq[TOVOLA_NU];
    TOVdata->Iso_r_arr[TOVdata->numpoints_actually_saved] =
        TOVola_eq[TOVOLA_R_ISO];
    TOVdata->numpoints_actually_saved++;

    /* Termination condition */
    if (TOVola_do_we_terminate(current_position, TOVola_eq, TOVdata)) {
      CCTK_VINFO("Finished Integration at position %.6e with Mass %.14e",
                 current_position, TOVola_eq[TOVOLA_MASS]);
      break;
    }
  }

  /* Cleanup */
  gsl_odeiv2_driver_free(driver);

  // Data in TOVdata->*_arr are stored at r=TOVdata->initial_ode_step_size > 0
  // up to the stellar surface. However, we may need data at r=0, which would
  // require extrapolation. To prevent that, we copy INTERP_BUFFER data points
  // from r>0 to r<0 so that we can always interpolate.
  CCTK_REAL *restrict tmp =
      malloc(sizeof(CCTK_REAL) *
             (TOVdata->numpoints_actually_saved + NEGATIVE_R_INTERP_BUFFER));

  extend_to_negative_r(TOVdata->rSchw_arr, -1.0, tmp, TOVdata);
  extend_to_negative_r(TOVdata->rho_energy_arr, +1.0, tmp, TOVdata);
  extend_to_negative_r(TOVdata->rho_baryon_arr, +1.0, tmp, TOVdata);
  extend_to_negative_r(TOVdata->P_arr, +1.0, tmp, TOVdata);
  extend_to_negative_r(TOVdata->M_arr, +1.0, tmp, TOVdata);
  extend_to_negative_r(TOVdata->nu_arr, +1.0, tmp, TOVdata);
  extend_to_negative_r(TOVdata->Iso_r_arr, -1.0, tmp, TOVdata);

  free(tmp);
  TOVdata->numpoints_actually_saved += NEGATIVE_R_INTERP_BUFFER;

  /* Allocate and populate TOVola_ID_persist_struct arrays */
  initialize_ID_persist_data(TOVola_ID_persist, TOVdata);

  TOVola_ID_persist->numpoints_arr = TOVdata->numpoints_actually_saved;

  /* Normalize and set data */
  TOVola_Normalize_and_set_data_integrated(
      TOVdata, TOVola_ID_persist->r_Schw_arr, TOVola_ID_persist->rho_energy_arr,
      TOVola_ID_persist->rho_baryon_arr, TOVola_ID_persist->P_arr,
      TOVola_ID_persist->M_arr, TOVola_ID_persist->expnu_arr,
      TOVola_ID_persist->exp4phi_arr, TOVola_ID_persist->r_iso_arr);

  /* Free raw data as it's no longer needed */
  free_tovola_data(TOVdata);

  /* Now to interp, and finalize the grid. */
  CCTK_REAL TOVola_Rbar =
      TOVola_ID_persist->r_iso_arr[TOVola_ID_persist->numpoints_arr - 1];
  CCTK_REAL TOVola_Mass =
      TOVola_ID_persist->M_arr[TOVola_ID_persist->numpoints_arr - 1];

  // Now for the actual grid placements. Go over all grid points
  CCTK_LOOP3_ALL(TOVola_Grid_Interp, cctkGH, i, j, k) {
    CCTK_INT i3d = CCTK_GFINDEX3D(cctkGH, i, j, k); // 3D index
    CCTK_INT i3d_vx =
        CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 0); // index for velocity_x
    CCTK_INT i3d_vy =
        CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 1); // index for velocity_y
    CCTK_INT i3d_vz =
        CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 2); // index for velocity_z
    CCTK_REAL TOVola_r_iso =
        sqrt((x[i3d] * x[i3d]) + (y[i3d] * y[i3d]) +
             (z[i3d] * z[i3d])); // magnitude of r on the grid
    CCTK_REAL TOVola_rho_energy, TOVola_rho_baryon, TOVola_P, TOVola_M,
        TOVola_expnu, TOVola_exp4phi; // Declare TOV quantities

    // Need to respect the initial parameters in Hydrobase and ADMbase. We will
    // do this by making local variables that we will write to the gridfunctions
    // after interpolation. We will start by assuming we are outside the star,
    // so we can reduce the number of times we interpolate the TOV data.
    CCTK_REAL TOVola_rSchw_outside =
        (TOVola_r_iso + TOVola_Mass) +
        TOVola_Mass * TOVola_Mass /
            (4.0 *
             TOVola_r_iso); // Need to know what rSchw is at our current point.
    CCTK_REAL TOVola_rho_local = 0.0;
    CCTK_REAL TOVola_press_local = 0.0;
    CCTK_REAL TOVola_eps_local = 0.0;
    CCTK_REAL TOVola_alp_local = sqrt(
        1 - 2 * TOVola_Mass / TOVola_rSchw_outside); // Goes to Schwarschild
    CCTK_REAL TOVola_gxx_local =
        pow((TOVola_rSchw_outside / TOVola_r_iso), 2.0);
    CCTK_REAL TOVola_gyy_local = TOVola_gxx_local;
    CCTK_REAL TOVola_gzz_local = TOVola_gxx_local;
    if (TOVola_r_iso < TOVola_Rbar) { // If we are INSIDE the star, we need to
                                      // interpollate the data to the grid.
      TOVola_TOV_interpolate_1D(
          TOVola_r_iso, TOVola_Interpolation_Stencil,
          TOVola_Max_Interpolation_Stencil, TOVola_ID_persist->numpoints_arr,
          TOVola_ID_persist->r_Schw_arr, TOVola_ID_persist->rho_energy_arr,
          TOVola_ID_persist->rho_baryon_arr, TOVola_ID_persist->P_arr,
          TOVola_ID_persist->M_arr, TOVola_ID_persist->expnu_arr,
          TOVola_ID_persist->exp4phi_arr, TOVola_ID_persist->r_iso_arr,
          &TOVola_rho_energy, &TOVola_rho_baryon, &TOVola_P, &TOVola_M,
          &TOVola_expnu, &TOVola_exp4phi);
      TOVola_rho_local = TOVola_rho_baryon;
      TOVola_press_local = TOVola_P;
      TOVola_eps_local = (TOVola_rho_energy / (TOVola_rho_baryon)) - 1.0;
      if (TOVola_eps_local < 0) {
        TOVola_eps_local = 0.0;
      }
      TOVola_alp_local = sqrt(TOVola_expnu); // This is the lapse
      TOVola_gxx_local = TOVola_exp4phi; // This is the values for the metric in
                                         // the coordinates we chose.
      TOVola_gyy_local = TOVola_gxx_local;
      TOVola_gzz_local = TOVola_gxx_local;
    }

    if (CCTK_EQUALS(initial_hydro, "TOVola")) {
      rho[i3d] = TOVola_rho_local;
      press[i3d] = TOVola_press_local;
      eps[i3d] = TOVola_eps_local;
      w_lorentz[i3d] = 1.0;
      vel[i3d_vx] = initial_velx;
      vel[i3d_vy] = initial_vely;
      vel[i3d_vz] = initial_velz;
    }
    if (CCTK_EQUALS(initial_data, "TOVola")) {
      gxx[i3d] = TOVola_gxx_local;
      gyy[i3d] = TOVola_gyy_local;
      gzz[i3d] = TOVola_gzz_local;
      gxy[i3d] = 0.0;
      gxz[i3d] = 0.0;
      gyz[i3d] = 0.0;
      kxx[i3d] = 0.0; // Curvature is zero for our slice
      kyy[i3d] = 0.0;
      kzz[i3d] = 0.0;
      kxy[i3d] = 0.0;
      kxz[i3d] = 0.0;
      kyz[i3d] = 0.0;
    }
    if (CCTK_EQUALS(initial_shift, "TOVola")) {
      betax[i3d] = initial_shiftx;
      betay[i3d] = initial_shifty;
      betaz[i3d] = initial_shiftz;
    }
    if (CCTK_EQUALS(initial_lapse, "TOVola")) {
      alp[i3d] = TOVola_alp_local;
    }
  }
  CCTK_ENDLOOP3_ALL(TOVola_Grid_Interp);

  // This is to populate time levels.
  // Luckily, this is a static solution, so the logic isn't too complicated.
  // Just copy data to the other timelevels
  switch (TOVola_TOV_Populate_Timelevels) {
  case 3:
    CCTK_LOOP3_ALL(TOVola_Pop2, cctkGH, i, j, k) {
      CCTK_INT i3d = CCTK_GFINDEX3D(cctkGH, i, j, k); // 3D index
      CCTK_INT i3d_vx =
          CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 0); // index for velocity_x
      CCTK_INT i3d_vy =
          CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 1); // index for velocity_y
      CCTK_INT i3d_vz =
          CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 2); // index for velocity_z
      // Copy Over
      if (CCTK_EQUALS(initial_hydro, "TOVola")) {
        rho_p_p[i3d] = rho[i3d];
        press_p_p[i3d] = press[i3d];
        eps_p_p[i3d] = eps[i3d];
        w_lorentz_p_p[i3d] = w_lorentz[i3d];
        vel_p_p[i3d_vx] = vel[i3d_vx];
        vel_p_p[i3d_vy] = vel[i3d_vy];
        vel_p_p[i3d_vz] = vel[i3d_vz];
      }
      if (CCTK_EQUALS(initial_data, "TOVola")) {
        gxx_p_p[i3d] = gxx[i3d];
        gyy_p_p[i3d] = gyy[i3d];
        gzz_p_p[i3d] = gzz[i3d];
        gxy_p_p[i3d] = gxy[i3d];
        gxz_p_p[i3d] = gxz[i3d];
        gyz_p_p[i3d] = gyz[i3d];
      }
    }
    CCTK_ENDLOOP3_ALL(TOVola_Pop2);
    // fall through
  case 2:
    CCTK_LOOP3_ALL(TOVola_Pop1, cctkGH, i, j, k) {
      CCTK_INT i3d = CCTK_GFINDEX3D(cctkGH, i, j, k); // 3D index
      CCTK_INT i3d_vx =
          CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 0); // index for velocity_x
      CCTK_INT i3d_vy =
          CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 1); // index for velocity_y
      CCTK_INT i3d_vz =
          CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 2); // index for velocity_z
      // Copy Over
      if (CCTK_EQUALS(initial_hydro, "TOVola")) {
        rho_p[i3d] = rho[i3d];
        press_p[i3d] = press[i3d];
        eps_p[i3d] = eps[i3d];
        w_lorentz_p[i3d] = w_lorentz[i3d];
        vel_p[i3d_vx] = vel[i3d_vx];
        vel_p[i3d_vy] = vel[i3d_vy];
        vel_p[i3d_vz] = vel[i3d_vz];
      }
      if (CCTK_EQUALS(initial_data, "TOVola")) {
        gxx_p[i3d] = gxx[i3d];
        gyy_p[i3d] = gyy[i3d];
        gzz_p[i3d] = gzz[i3d];
        gxy_p[i3d] = gxy[i3d];
        gxz_p[i3d] = gxz[i3d];
        gyz_p[i3d] = gyz[i3d];
      }
    }
    CCTK_ENDLOOP3_ALL(TOVola_Pop1);
    // fall through
  case 1:
    break;
  default:
    CCTK_ERROR("Unsupported number of TOVola_TOV_Populate_Timelevels");
    break;
  }

  free_ID_persist_data(TOVola_ID_persist);
}
