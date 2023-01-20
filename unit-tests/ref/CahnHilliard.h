#include "mechanoChemFEM.h"
#include <math.h> /* exp */
#include <cmath>  /* sqrt */
template <int dim>
class CahnHilliard : public mechanoChemFEM<dim>
{
public:
  CahnHilliard();
  // this is a overloaded function
  void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim> &fe_values, Table<1, Sacado::Fad::DFad<double>> &R, Table<1, Sacado::Fad::DFad<double>> &ULocal, Table<1, double> &ULocalConv);
  ParameterHandler *params;
  void apply_initial_condition();
  void solve_ibvp();
  ConstraintMatrix *constraints;

  int iter_count = 0;
  bool increase_ts = false;

  int c_dof = 0, mu_dof = 1, phi_dof = 2, zeta_dof = 3; // vacancy, mu, Bi layer (constraint), not used (constraint)
};
template <int dim>
CahnHilliard<dim>::CahnHilliard()
{
  constraints = this->constraints_mechanoChemFEM;
  // This let you use one params to get all parameters pre-defined in the mechanoChemFEM
  params = this->params_mechanoChemFEM;
  params->enter_subsection("Physics");
  params->declare_entry("omega", "0", Patterns::Double());
  params->declare_entry("kappa", "0", Patterns::Double());
  params->declare_entry("M", "0", Patterns::Double());
  // init
  params->declare_entry("perturb", "0", Patterns::Double());
  params->declare_entry("mean", "0", Patterns::Double());
  // mu
  params->declare_entry("mu_fun", "0", Patterns::Integer());
  params->declare_entry("constant_M","false",Patterns::Bool());
  params->declare_entry("increase_time","false",Patterns::Bool());
  params->declare_entry("ep", "0", Patterns::Double());
  params->declare_entry("tt", "0", Patterns::Double());
  params->declare_entry("order", "1", Patterns::Integer());
  // flux
  params->declare_entry("flux", "0.1", Patterns::Double());
  params->declare_entry("flux_width", "0.1", Patterns::Double());
  // mobility
  params->declare_entry("constant_M","false",Patterns::Bool());
  params->declare_entry("aa", "0.033", Patterns::Double());  // eV
  params->declare_entry("bb", "0.344", Patterns::Double());  // eV
  params->declare_entry("cc", "-0.377", Patterns::Double()); // eV
  params->leave_subsection();

  // Declare the parameters before load it
  this->load_parameters("parameters.prm");

  // define main fields from parameter file.
  this->define_primary_fields();
  // Set up the ibvp.
  this->init_ibvp();
}

long factorial(const int n)
{
  long f = 1;
  for (int i=1; i<=n; ++i)
      f *= i;
  return f;
}

template <int dim>
void CahnHilliard<dim>::apply_initial_condition()
{
  std::cout << "=========== applying initial condition (new) ===========\n";
  params->enter_subsection("Geometry");
  double x_min = params->get_double("x_min");
  double x_max = params->get_double("x_max");
  double y_min = params->get_double("y_min");
  double y_max = params->get_double("y_max");
  params->leave_subsection();

  params->enter_subsection("Physics");
  double perturb = params->get_double("perturb");
  double mean = params->get_double("mean");
  params->leave_subsection();

  // 0, 1, 2, 3 for different dof
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
  {
    if (cell->subdomain_id() == this->this_mpi_process)
    {
      Point<dim> center = cell->center();
      hp::FEValues<dim> hp_fe_values(this->fe_collection, this->q_collection, update_values | update_quadrature_points);
      hp_fe_values.reinit(cell);
      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
      const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      std::vector<unsigned int> local_dof_indices(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);

      int vertex_id = 0;
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        int ck = fe_values.get_fe().system_to_component_index(i).first;
        if (ck == c_dof)
        {
          this->solution_prev(local_dof_indices[i]) =  mean + perturb * (static_cast<double>(rand()) / (static_cast<double>(RAND_MAX)) - 0.5) * 2;
        }
        if (ck == mu_dof)
        {
          this->solution_prev(local_dof_indices[i]) = 0.0;
        }
        if (ck == phi_dof)
        {
          this->solution_prev(local_dof_indices[i]) = 0.0;
        }
        if (ck == zeta_dof)
        {
          this->solution_prev(local_dof_indices[i]) = 0.0;
        }
        if (ck == zeta_dof)
          vertex_id += 1;
      } // dofs_per_cell
    }   // this_mpi_process
  }     // cell

  this->solution_prev.compress(VectorOperation::insert); // potentially for parallel computing purpose.
  this->solution = this->solution_prev;

  constraints->clear();
  DoFTools::make_hanging_node_constraints(this->dof_handler, *constraints);


  //-------------------
  // add constraints and dirchle BC
  //------------------- 
  {
    hp::FEValues<dim> hp_fe_values(this->fe_collection, this->q_collection, update_values | update_quadrature_points);
    typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->subdomain_id() == this->this_mpi_process)
      { // parallel computing.
        // int cell_id = cell->active_cell_index();
        hp_fe_values.reinit(cell);
        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        std::vector<unsigned int> local_dof_indices(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
          // std::cout << cell_id << " ck " << ck << std::endl;
          if (ck == phi_dof or ck == zeta_dof)
          {
            auto globalDOF = local_dof_indices[i];
            constraints->add_line(globalDOF);
            constraints->set_inhomogeneity(globalDOF, 0.0);
          }
        }
      }
    }
  }

  constraints->close();
  std::cout << "===========  end of constraint =========== " << std::endl;
}

template <int dim>
void CahnHilliard<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim> &fe_values, Table<1, Sacado::Fad::DFad<double>> &R, Table<1, Sacado::Fad::DFad<double>> &ULocal, Table<1, double> &ULocalConv)
{

  params->enter_subsection("Geometry");
  double x_min = params->get_double("x_min");
  double x_max = params->get_double("x_max");
  double y_min = params->get_double("y_min");
  double y_max = params->get_double("y_max");
  // int num_elem_x = params->get_integer("num_elem_x");
  // int num_elem_y = params->get_integer("num_elem_y");
  params->leave_subsection();

  // evaluate primary fields
  params->enter_subsection("Physics");
  double M = params->get_double("M");
  double omega = params->get_double("omega");
  double kappa = params->get_double("kappa");
  // mu
  int mu_fun = params->get_integer("mu_fun");
  double ep = params->get_double("ep");
  double tt = params->get_double("tt");
  int order = params->get_integer("order");
  // flux
  double flux_width = params->get_double("flux_width");
  double flux = params->get_double("flux");
  // mobility
  bool constant_M = params->get_bool("constant_M");
  double aa = params->get_double("aa");
  double bb = params->get_double("bb");
  double cc = params->get_double("cc");
  params->leave_subsection();
  unsigned int n_q_points = fe_values.n_quadrature_points;

  dealii::Table<1, double> c_1_conv(n_q_points);
  dealii::Table<1, Sacado::Fad::DFad<double>> c_1(n_q_points), mu(n_q_points), phi(n_q_points), zeta(n_q_points);
  dealii::Table<2, Sacado::Fad::DFad<double>> c_1_grad(n_q_points, dim), mu_grad(n_q_points, dim), phi_grad(n_q_points, dim), zeta_grad(n_q_points, dim);
  dealii::Table<2, double> c_1_grad_conv(n_q_points, dim);

  evaluateScalarFunction<double, dim>(fe_values, c_dof, ULocalConv, c_1_conv);
  evaluateScalarFunctionGradient<double, dim>(fe_values, c_dof, ULocalConv, c_1_grad_conv);

  evaluateScalarFunction<Sacado::Fad::DFad<double>, dim>(fe_values, c_dof, ULocal, c_1);
  evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, dim>(fe_values, c_dof, ULocal, c_1_grad);

  evaluateScalarFunction<Sacado::Fad::DFad<double>, dim>(fe_values, mu_dof, ULocal, mu);
  evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, dim>(fe_values, mu_dof, ULocal, mu_grad);

  evaluateScalarFunction<Sacado::Fad::DFad<double>, dim>(fe_values, phi_dof, ULocal, phi);
  evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, dim>(fe_values, phi_dof, ULocal, phi_grad);

  evaluateScalarFunction<Sacado::Fad::DFad<double>, dim>(fe_values, zeta_dof, ULocal, zeta);
  evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, dim>(fe_values, zeta_dof, ULocal, zeta_grad);

  // evaluate diffusion and reaction term
  dealii::Table<1, Sacado::Fad::DFad<double>> rhs_mu(n_q_points);
  dealii::Table<2, Sacado::Fad::DFad<double>> j_c_1(n_q_points, dim), kappa_c_1_grad(n_q_points, dim);

  // scalar M
  // j_c_1 = table_scaling<dim, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double>>(mu_grad, -M); 

  double T = 300;           // K
  
  double a = -1.51629e+10 / 1e27 * 6.242e18;  // ev/nm3
  double b = 1.95270e+10  / 1e27 * 6.242e18;  // ev/nm3
  double c = -4.36414e+09 / 1e27 * 6.242e18;  // ev/nm3
  
  // double k = 1.380649e-23 * 6.242e18;        // m2 kg s-2 K-1 = J/K => eV/K;
  double k = 8.617333262e-5;                  // eV/K
  double V = 20.1207e-30  * 1e27;              // m^3

  // double aa = 5.226199007e-21 * 6.242e18;  // eV
  // double bb = 5.51601893e-20  * 6.242e18;  // eV
  // double cc = -6.03863884e-20 * 6.242e18;  // eV

  double nus = 1e13;
  // double DE =  0.0; //110.535 * 1.60218E-22; // mev - > J;
  
  Sacado::Fad::DFad<double> dGdc;
  Sacado::Fad::DFad<double> L = M;
  double d;
  double P;
  
  Sacado::Fad::DFad<double> Y;
  double ep0 = tt +    ep  * (1-2*tt);
  double ep1 = tt + (1-ep) * (1-2*tt);
  for (unsigned int q = 0; q < n_q_points; q++)
  { 
    // ------------------------------------------------------------------------------------------------
    if (mu_fun == 0)
    {
      this->pcout << "right" << "\n";
      if ((ep0 < c_1[q]) and (c_1[q] < ep1))
      {
        Y = (c_1[q] - tt) / (1 - 2*tt);
        dGdc = a + b + 2*c*Y +(k*T/V)*(log(Y) - log(1 - Y));   // J/m^3
      }
      else
      {
        if (c_1[q] <= ep0)
        {
          P = ep; // =  Y(ep0)
          this->pcout << "left" << "\n";
        }
        else if (ep1 <= c_1[q])
        {
          P = 1 - ep; // =  Y(ep1)
          this->pcout << "right" << "\n";
        }

        dGdc = a + b + 2*c*P + (k*T/V)*(log(P) - log(1 - P)); 
        if (order>=1)
        {
          dGdc += 2*c/(1-2*tt)*(c_1[q]-P);
          for (unsigned int l = 1; l <= (order-1); l++)
          {                                                                                            
            dGdc += (k*T/V)/(1-2*tt)*std::pow(-1,l-1)*factorial(l-1)*(std::pow(P,-l)-std::pow(1-P,-l)) * std::pow(c_1[q] - P,l)/factorial(l);
          }
        }

      }
      dGdc = dGdc * V;
    }  
    else if (mu_fun == 1)
    {
      d = 5e2;
      dGdc =  0.865 * (a + b + 2 * c * c_1[q] - (T * k * (d*exp(-d*c_1[q]) - d*exp(d*(c_1[q] - 1))))/V ) * V;  // J/m^3
    }

    if (!constant_M)
    {
      L = M * c_1[q] * (1 - c_1[q]) * std::exp(-(aa + bb * c_1[q] + cc * c_1[q] * c_1[q]) / (k * T)); // no unit
      // if (c_1[q] <= ep)
      // {
      //   L = -1e9;
      // }
      // else if ((1 - ep) <= c_1[q])
      // {
      //   L = -1e9;
      // }
    }
    // ------------------------------------------------------------------------------------------------

    rhs_mu[q] = - mu[q] + omega * dGdc;
    
    for (unsigned int j = 0; j < dim; j++)
    {
      kappa_c_1_grad(q, j) =  kappa * c_1_grad(q, j);
      j_c_1(q, j) = - L * mu_grad(q, j);
    }
  }

  // call residual functions
  this->ResidualEq.residualForDiffusionEq(fe_values, c_dof, R, c_1, c_1_conv, j_c_1);
  this->ResidualEq.residualForPoissonEq(fe_values, mu_dof, R, kappa_c_1_grad, rhs_mu);

  //-------------------
  // BC
  //-------------------
  double tol = 0.001;
  for (unsigned int faceID = 0; faceID < 2 * dim; faceID++)
  {
    if (cell->face(faceID)->at_boundary() == true)
    {
      if ((cell->face(faceID)->center()[1] <= (y_min + tol)) and
          (cell->face(faceID)->center()[0] >= ((x_max + x_min) / 2 - flux_width / 2)) and
          (cell->face(faceID)->center()[0] <= ((x_max + x_min) / 2 + flux_width / 2))) // bottom boundary 
      {
        flux = - flux;
      }
      else
      {
        flux = 0.0;
      }
    FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
    fe_face_values.reinit(cell, faceID);
    this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, c_dof, R, flux); // minus flux for inbound // plus flux for outbound
    }
  }

  const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
  for (unsigned int i0 = 0; i0 < dofs_per_cell; ++i0)
  {
    if (R[i0] != R[i0])
    {
      std::cout << " R[i0] " << R[i0] << std::endl;
      exit(0);
    }
  }

}

/**************************************
 * Implement adaptive time stepping
 *************************************/
template <int dim>
void CahnHilliard<dim>::solve_ibvp()
{

  params->enter_subsection("Physics");
  bool increase_time = params->get_bool("increase_time");
  params->leave_subsection();

  if (increase_ts and increase_time)
  {
    increase_ts = false;
    this->pcout << "Increasing timestep \n";
    this->current_dt *= std::pow(10., 0.25); // increase dt
  }

  int converged;
  while (true)
  {
    converged = this->nonlinearSolve(this->solution);
    if ((converged > -1))
    {
      break;
    }
    else
    {
      iter_count = 0;
      this->pcout << "Not converged or out of bounds, reduce dt." << std::endl;
      this->current_dt /= std::pow(10., 0.25);
      this->solution = this->solution_prev;
    }
  }
  if (converged < 4)
  {
    iter_count++;
    if (converged < 3)
    {
      iter_count++;
      iter_count++; // Increase time step sooner for faster convergence
    }
  }
  else
  {
    iter_count = 0;
  }
  // Check if the current_time is a multiple of the doubled timestep
  // bool multiple = (std::fmod(this->current_time,2.*this->current_dt) < 1.e-8);
  if (iter_count > 2 && this->current_dt < 0.1)
  {
    increase_ts = true;
    iter_count = 0;
  }

  this->solution_prev = this->solution;
}
