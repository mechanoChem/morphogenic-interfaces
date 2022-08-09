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

  int c_dof = 0, mu_dof = 1, phi_dof = 2, zeta_dof = 3; // c, mu, LLZO (constraint), Bi (constraint)
};
template <int dim>
CahnHilliard<dim>::CahnHilliard()
{
  constraints = this->constraints_mechanoChemFEM;
  // This let you use one params to get all parameters pre-defined in the mechanoChemFEM
  params = this->params_mechanoChemFEM;
  params->enter_subsection("Concentration");
  params->declare_entry("c_ini", "0", Patterns::Double());
  params->declare_entry("omega", "0", Patterns::Double());
  params->declare_entry("c_alpha", "0", Patterns::Double());
  params->declare_entry("c_beta", "0", Patterns::Double());
  params->declare_entry("kappa", "0", Patterns::Double());
  params->declare_entry("M", "0", Patterns::Double());
  params->declare_entry("phi_0", "0", Patterns::Double());
  params->declare_entry("zeta_0", "0", Patterns::Double());
  // boundary flux
  params->declare_entry("max_flux", "0.1", Patterns::Double());
  params->declare_entry("min_flux", "0.01", Patterns::Double());
  params->declare_entry("spike_loc", "0.5", Patterns::Double());
  params->declare_entry("sharpness", "50.0", Patterns::Double());
  params->leave_subsection();

  // Declear the parameters before load it
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
        // if (ck==c_dof) this->solution_prev(local_dof_indices[i]) = exp(-10.0 * cell->vertex(vertex_id)[0]) * exp(-10.0 * (1.0-cell->vertex(vertex_id)[1]));
        //  if (ck==c_dof) this->solution_prev(local_dof_indices[i]) = exp(-5.0 * cell->vertex(vertex_id)[0]) * exp(-5.0 * (1.0-cell->vertex(vertex_id)[1]));
        if (ck == c_dof)
          this->solution_prev(local_dof_indices[i]) = 0.01; //exp(-5.0 * cell->vertex(vertex_id)[0]) * exp(-5.0 * cell->vertex(vertex_id)[1]);
        // if (ck == c_dof and this->solution_prev(local_dof_indices[i]) < 0.01)
          // this->solution_prev(local_dof_indices[i]) = 0.01; // minimum limit
        if (ck == mu_dof)
          this->solution_prev(local_dof_indices[i]) = 0.0;
        if (ck == phi_dof)
          this->solution_prev(local_dof_indices[i]) = exp(-30.0 * cell->vertex(vertex_id)[1]); // phi : LLZO
        if (ck == zeta_dof)
          this->solution_prev(local_dof_indices[i]) = exp(-30.0 * (y_max - cell->vertex(vertex_id)[1])); // zeta: Bi
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
  // add constraints
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
  double c_alpha = params->get_double("c_alpha");
  double c_beta = params->get_double("c_beta");
  double kappa = params->get_double("kappa");
  double phi_0 = params->get_double("phi_0");
  double zeta_0 = params->get_double("zeta_0");
  double max_flux = params->get_double("max_flux");
  double min_flux = params->get_double("min_flux");
  double spike_loc = params->get_double("spike_loc");
  double sharpness = params->get_double("sharpness");
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

  kappa_c_1_grad = table_scaling<dim, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double>>(c_1_grad, kappa);

  // scalar M
   j_c_1=table_scaling<dim, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double> >(mu_grad,-M);//-D_1*c_1_grad

  for (unsigned int q = 0; q < n_q_points; q++)
    rhs_mu[q] = (omega * (log(c_1[q] / (1.0 - c_1[q]))) + 2.5 * omega * (1.0 - 2.0 * c_1[q]))- mu[q];
        // omega=A, B=2.5*A
  // call residual functions
  this->ResidualEq.residualForDiffusionEq(fe_values, c_dof, R, c_1, c_1_conv, j_c_1);
  this->ResidualEq.residualForPoissonEq(fe_values, mu_dof, R, kappa_c_1_grad, rhs_mu);

  //-------------------
  // BC
  //-------------------
  double range = sqrt((max_flux - min_flux) / sharpness);
  // std::cout << " max_flux = " << max_flux << std::endl;
  // std::cout << " min_flux = " << min_flux << std::endl;
  // std::cout << " sharpness = " << sharpness << std::endl;
  // std::cout << " rang = " << range << std::endl;
  double flux = min_flux;
  for (unsigned int faceID = 0; faceID < 2 * dim; faceID++)
  {
    if (cell->face(faceID)->at_boundary() == true)
    {
      // if(cell->face(faceID)->center()[0] <= 0.001 and cell->face(faceID)->center()[1] >= 0.97)
      if (cell->face(faceID)->center()[1] <= (y_min + 0.001))
      {
        if (((spike_loc - range ) <= cell->face(faceID)->center()[0]) and (cell->face(faceID)->center()[0] <= (spike_loc + range )))
        {
          flux = -sharpness * pow(cell->face(faceID)->center()[0] - spike_loc, 2) + max_flux;
        }
        else
        {
          flux = min_flux;
        }
        // std::cout << " c_1[q] " << c_1[0].val() << " " << c_1[1].val() << " " << c_1[2].val() << " " << c_1[3].val() << std::endl;
          FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
          fe_face_values.reinit(cell, faceID);
          this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, c_dof, R, -flux);
          // std::cout << "flux = " << -flux << std::endl;
      }
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

template <int dim>
void InitialConditions<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const
{
  Assert(values.size() == 2, ExcDimensionMismatch(values.size(), 2));
  values(1) = 0;
  values(0) = 0.5 + 0.04 * (static_cast<double>(rand()) / (static_cast<double>(RAND_MAX)) - 0.5);
}
template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;
