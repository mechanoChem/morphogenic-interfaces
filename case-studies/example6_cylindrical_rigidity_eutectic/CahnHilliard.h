/*
zhenlin wang 2019
*CahnHilliard
*/
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
  ConstraintMatrix *constraints;

  int c_dof = 0, mu_dof = 1, phi_dof = 2, zeta_dof = 3; // c, mu, Bi layer (constraint), not used (constraint)
};
template <int dim>
CahnHilliard<dim>::CahnHilliard()
{
  constraints = this->constraints_mechanoChemFEM;
  // This let you use one params to get all parameters pre-defined in the mechanoChemFEM
  params = this->params_mechanoChemFEM;
  params->enter_subsection("Concentration");
  params->declare_entry("omega", "0", Patterns::Double());
  params->declare_entry("kappa", "0", Patterns::Double());
  params->declare_entry("M", "0", Patterns::Double());
  params->declare_entry("c_ini", "0.0", Patterns::Double());
  // gap
  params->declare_entry("ww", "0.25", Patterns::Double());
  params->declare_entry("hh", "0.25", Patterns::Double());
  params->declare_entry("dd", "0.1", Patterns::Double());
  params->declare_entry("rate", "20.0", Patterns::Double());
  params->declare_entry("phi_amp", "0.0", Patterns::Double());
  params->declare_entry("zeta_amp", "0.0", Patterns::Double());
  // boundary flux
  params->declare_entry("max_flux", "0.1", Patterns::Double());
  params->leave_subsection();

  // Declare the parameters before load it
  this->load_parameters("parameters.prm");

  // define main fields from parameter file.
  this->define_primary_fields();
  // Set up the ibvp.
  this->init_ibvp();
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

  params->enter_subsection("Concentration");
  double ww = params->get_double("ww");
  double hh = params->get_double("hh");
  double dd = params->get_double("dd");
  double c_ini = params->get_double("c_ini");
  double rate = params->get_double("rate");
  double phi_amp = params->get_double("phi_amp");
  double zeta_amp = params->get_double("zeta_amp");
  params->leave_subsection();

  double a1,b1, a2,b2, a3,b3;
  double D1, D2, D3;
  double C1, C2, C3;

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
        a1 = 1. * hh;
        a2 = 5. * hh; 
        a3 = 9. * hh;
        b1 = -1. * ww;
        b2 = -5. * ww;
        b3 = -9. * ww;
        D1 = abs(a1 * cell->vertex(vertex_id)[0] + b1 * cell->vertex(vertex_id)[1] + a1*b1) / sqrt(pow(a1,2) + pow(b1,2));
        D2 = abs(a2 * cell->vertex(vertex_id)[0] + b2 * cell->vertex(vertex_id)[1] + a2*b2) / sqrt(pow(a2,2) + pow(b2,2));
        D3 = abs(a3 * cell->vertex(vertex_id)[0] + b3 * cell->vertex(vertex_id)[1] + a3*b3) / sqrt(pow(a3,2) + pow(b3,2));
        C1 = 1.0 - 1.0 / (1.0 + exp((rate / dd) * (-D1 + dd)));
        C2 = 1.0 - 1.0 / (1.0 + exp((rate / dd) * (-D2 + dd)));
        C3 = 1.0 - 1.0 / (1.0 + exp((rate / dd) * (-D3 + dd)));
        if (ck == c_dof)
          this->solution_prev(local_dof_indices[i]) = c_ini; // exp(-5.0 * cell->vertex(vertex_id)[0]) * exp(-5.0 * cell->vertex(vertex_id)[1]);
        if (ck == mu_dof)
          this->solution_prev(local_dof_indices[i]) = 0.0;
        if (ck == phi_dof)
        {
          this->solution_prev(local_dof_indices[i]) = (C1 + C2 + C3) * phi_amp;
        }
        if (ck == zeta_dof)
        {
          this->solution_prev(local_dof_indices[i]) = (C1 + C2 + C3) * zeta_amp;
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
  // add constraints and BC
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
  params->enter_subsection("Concentration");
  double M = params->get_double("M");
  double omega = params->get_double("omega");
  double kappa = params->get_double("kappa");
  double max_flux = params->get_double("max_flux");
  double ww = params->get_double("ww");
  double hh = params->get_double("hh");
  double dd = params->get_double("dd");
  params->leave_subsection();
  unsigned int n_q_points = fe_values.n_quadrature_points;

  double T = 300; // K
  double a0 = -1.51629e+10 / 1e27 * 6.242e18; // ev/nm3
  double b0 = 1.95270e+10 / 1e27 * 6.242e18;  // ev/nm3
  double c0 = -4.36414e+09 / 1e27 * 6.242e18; // ev/nm3

  // double k = 1.380649e-23 * 6.242e18;        // m2 kg s-2 K-1 = J/K => eV/K;
  double k = 8.617333262e-5;     // eV/K
  double V = 0.0201207;          // nm^3   = 20.1207e-30 * 1e27;
  double d = 50;



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
  
  Sacado::Fad::DFad<double> f_c;

  // scalar M
  // j_c_1 = table_scaling<dim, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double>>(mu_grad, -M); 

  for (unsigned int q = 0; q < n_q_points; q++)
  {
    // rhs_mu[q] = -mu[q] + (1.0 - phi[q]) * omega * (4 * c_1[q] * c_1[q] * c_1[q] - 6 * c_1[q] * c_1[q] + 2 * c_1[q] - 57.5646 * std::pow(10, -50.0 * c_1[q]) + 5.75646 * std::pow(10, -49) * std::pow(10, 50.0 * c_1[q]));
    
    f_c = 0.865 * (a0 + b0 + 2 * c0 * c_1[q] + (k * T/ V * d) * (exp(d * (c_1[q] - 1)) - exp(-d * c_1[q]))) ;   // eV/nm^3

    // Cahn-Hilliard    
    rhs_mu[q] = -mu[q] + (1.0 - phi[q]) * f_c;   

    for (unsigned int j = 0; j < dim; j++)
    {
      kappa_c_1_grad(q, j) = (1.0 - phi[q]) * kappa * c_1_grad(q, j);
      j_c_1(q, j) = -M * (1.0 - zeta[q]) * mu_grad(q, j);
    }
  }

  // call residual functions
  this->ResidualEq.residualForDiffusionEq(fe_values, c_dof, R, c_1, c_1_conv, j_c_1);
  this->ResidualEq.residualForPoissonEq(fe_values, mu_dof, R, kappa_c_1_grad, rhs_mu);

  //-------------------
  // BC
  //-------------------
  double flux;
  double tol = 0.001;
  for (unsigned int faceID = 0; faceID < 2 * dim; faceID++)
  {
    if (cell->face(faceID)->at_boundary() == true)
    {
      flux = 0.0;

      if (cell->face(faceID)->center()[1] <= (y_min +tol)) // bottom boundary 
      {
        if ((cell->face(faceID)->center()[0] >= (1 * ww + dd)) and
            (cell->face(faceID)->center()[0] <= (5 * ww - dd)))  //in gap
        {
          flux = - max_flux;
        }
        else if ((cell->face(faceID)->center()[0] >= (5 * ww + dd)) and
                 (cell->face(faceID)->center()[0] <= (9 * ww - dd))) //in gap
        {
          flux = - max_flux;
        }
      }
      // else if ((cell->face(faceID)->center()[0] >= (x_max - tol)) and
      //          (cell->face(faceID)->center()[1] >= (r2 + tol))) // right boundary above Bi       
      // {
      //   flux = max_flux / 2;
      // }

      // if ()

      // if (cell->face(faceID)->center()[1] <= (y_min +tol))
      //   flux = - max_flux;

      FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
      fe_face_values.reinit(cell, faceID);
      this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, c_dof, R, flux); // minus flux for inbound // plus fulx for outbound
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

