/*************************************
 * Cahn-Hilliard, two fields c1, c2
 * Mostafa Shojaei
 * Sun May 21 12:17:13 PDT 2023
 *************************************/

#include "mechanoChemFEM.h"
#include "extraFunctions.h"
// #include "supplementary/supplementaryFunctions.h" //table scaling ?
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

  int c1_dof = 0, mu1_dof = 1, c2_dof = 2, mu2_dof = 3; 
};

template <int dim>
CahnHilliard<dim>::CahnHilliard()
{
  constraints = this->constraints_mechanoChemFEM;
  // This let you use one params to get all parameters pre-defined in the mechanoChemFEM
  params = this->params_mechanoChemFEM;
  params->enter_subsection("Concentration");
  params->declare_entry("c1_ini", "0.25", Patterns::Double());
  params->declare_entry("c2_ini", "0.25", Patterns::Double());
  params->declare_entry("c1_per", "0.01", Patterns::Double());
  params->declare_entry("c2_per", "0.01", Patterns::Double());
  //
  params->declare_entry("mobility_1","0",Patterns::Double() );
  params->declare_entry("mobility_2","0",Patterns::Double() );
  params->declare_entry("kappa_1","0",Patterns::Double() );
  params->declare_entry("kappa_2","0",Patterns::Double() );
  // 
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
  double c1_ini = params->get_double("c1_ini");
  double c2_ini = params->get_double("c2_ini");
  double c1_per = params->get_double("c1_per");
  double c2_per = params->get_double("c2_per");
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

      int vertex_id = -1;
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        int ck = fe_values.get_fe().system_to_component_index(i).first;
        if (ck == 0) vertex_id += 1;
        if (ck == c1_dof)
          this->solution_prev(local_dof_indices[i]) = c1_ini + c1_per*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
        if (ck == mu1_dof)
          this->solution_prev(local_dof_indices[i]) = 0.0;
        if (ck == c2_dof)
          this->solution_prev(local_dof_indices[i]) = c2_ini + c2_per*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
        if (ck == mu2_dof)
          this->solution_prev(local_dof_indices[i]) = 0.0;

      } // dofs_per_cell
    }   // this_mpi_process
  }     // cell

  this->solution_prev.compress(VectorOperation::insert); // potentially for parallel computing purpose.
  this->solution = this->solution_prev;

  constraints->clear();

  // DoFTools::make_hanging_node_constraints(this->dof_handler, *constraints);
  // //-------------------
  // // add constraints and BC
  // //-------------------
  // {
  //   hp::FEValues<dim> hp_fe_values(this->fe_collection, this->q_collection, update_values | update_quadrature_points);
  //   typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
  //   for (; cell != endc; ++cell)
  //   {
  //     if (cell->subdomain_id() == this->this_mpi_process)
  //     { // parallel computing.
  //       // int cell_id = cell->active_cell_index();
  //       hp_fe_values.reinit(cell);
  //       const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

  //       const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
  //       std::vector<unsigned int> local_dof_indices(dofs_per_cell);
  //       cell->get_dof_indices(local_dof_indices);
  //       for (unsigned int i = 0; i < dofs_per_cell; ++i)
  //       {
  //         const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
  //         // std::cout << cell_id << " ck " << ck << std::endl;
  //         // if (ck == c1_dof)
  //         // {
  //         //   auto globalDOF = local_dof_indices[i];
  //         //   constraints->add_line(globalDOF);
  //         //   constraints->set_inhomogeneity(globalDOF, 0.0);
  //         // }
  //       }
  //     }
  //   }
  // }

  constraints->close();
  // std::cout << "===========  end of constraint =========== " << std::endl;
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
  double M1 = params->get_double("mobility_1");
	double M2 = params->get_double("mobility_2");
	double k1 = params->get_double("kappa_1");
	double k2 = params->get_double("kappa_2");
  params->leave_subsection();

  std::vector<double> x = {0.1, 0.90, 0.10, 0.45};
  std::vector<double> y = {0.10, 0.10, 0.90, 0.45};
  std::vector<double> v = {-1.0, -1.0, -1.0, 0.2};
  std::vector<double> r = {8.0, 8.0, 8.0, 1.0};
  double a = 1.0;
  // double m = -1.0; //= min v_val

  // --------------------------------------------------------------------

 	unsigned int n_q_points= fe_values.n_quadrature_points;

  //define fields
	dealii::Table<1,double>  c1_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > c1(n_q_points), mu1(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> >  c1_grad(n_q_points, dim), mu1_grad(n_q_points, dim);
	
	dealii::Table<1,double>  c2_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > c2(n_q_points), mu2(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> >  c2_grad(n_q_points, dim), mu2_grad(n_q_points, dim);
  //evaluate fields
	evaluateScalarFunction<double,dim>(fe_values, c1_dof, ULocalConv, c1_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c1_dof, ULocal, c1);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c1_dof, ULocal, c1_grad);
	evaluateScalarFunction<double,dim>(fe_values, c2_dof, ULocalConv, c2_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c2_dof, ULocal, c2);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c2_dof, ULocal, c2_grad);
	
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, mu1_dof, ULocal, mu1);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, mu1_dof, ULocal, mu1_grad);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, mu2_dof, ULocal, mu2);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, mu2_dof, ULocal, mu2_grad);
	
	
	//evaluate diffusion and reaction term
	dealii::Table<1,Sacado::Fad::DFad<double> > rhs_mu1(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_c1(n_q_points, dim), kappa_c1_grad(n_q_points, dim);
	dealii::Table<1,Sacado::Fad::DFad<double> > rhs_mu2(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_c2(n_q_points, dim), kappa_c2_grad(n_q_points, dim);


  Sacado::Fad::DFad<double> F_c1, F_c2;

  double d = .4;
  double s = .7;
	for(unsigned int q=0; q < n_q_points; q++){

    copy_to_q(j_c1, -M1 * copy_from_q(mu1_grad, q), q);  //j_c1[q][i] = -M_1 * c_1_grad[q][i], for i=1:dim
    copy_to_q(j_c2, -M2 * copy_from_q(mu2_grad, q), q);
    copy_to_q(kappa_c1_grad, k1 * copy_from_q(c1_grad, q), q);
    copy_to_q(kappa_c2_grad, k2 * copy_from_q(c2_grad, q), q);
		
    F_c1 = 0.0;
    F_c2 = 0.0;
    for (unsigned int i = 0; i < 4; i++){
      F_c1 += -r[i] * v[i] * exp(-r[i]*(std::pow(c1[q] - x[i],2) + std::pow(c2[q] - y[i],2))) * (2*c1[q] - 2*x[i]);
      F_c2 += -r[i] * v[i] * exp(-r[i]*(std::pow(c1[q] - x[i],2) + std::pow(c2[q] - y[i],2))) * (2*c2[q] - 2*y[i]);
    }
    F_c1 = a * F_c1;
		F_c2 = a * F_c2;

    // F_c1 = 6*d/std::pow(s,4)*(c1[q]*c1[q]+c2[q]*c2[q])*c1[q]-6*d/std::pow(s,3)*c1[q]*c2[q]-3*d/std::pow(s,2)*c1[q];
    // F_c2 = 6*d/std::pow(s,4)*(c1[q]*c1[q]+c2[q]*c2[q])*c2[q]+3*d/std::pow(s,3)*c2[q]*c2[q]-3*d/std::pow(s,2)*c2[q];
		
    rhs_mu1[q] = F_c1 - mu1[q];
    rhs_mu2[q] = F_c2 - mu2[q];

	}
	
	//call residual functions
	this->ResidualEq.residualForDiffusionEq(fe_values, c1_dof, R, c1, c1_conv, j_c1);
  this->ResidualEq.residualForDiffusionEq(fe_values, c2_dof, R, c2, c2_conv, j_c2);

	this->ResidualEq.residualForPoissonEq(fe_values, mu1_dof, R, kappa_c1_grad, rhs_mu1);
	this->ResidualEq.residualForPoissonEq(fe_values, mu2_dof, R, kappa_c2_grad, rhs_mu2);
	
  //-------------------
  // BC
  //-------------------
  // for (unsigned int faceID = 0; faceID < 2 * dim; faceID++)
  // {
  //   if (cell->face(faceID)->at_boundary() == true)
  //   {
  //     if (cell->face(faceID)->center()[1] <= (y_min +tol)) // bottom boundary 
  //     {
  //       if ((cell->face(faceID)->center()[0] >= (1 * ww + dd)) and
  //           (cell->face(faceID)->center()[0] <= (5 * ww - dd)))  //in gap
  //       {
  //         flux = - max_flux;
  //       }
  //       else if ((cell->face(faceID)->center()[0] >= (5 * ww + dd)) and
  //                (cell->face(faceID)->center()[0] <= (9 * ww - dd))) //in gap
  //       {
  //         flux = - max_flux;
  //       }
  //     }
  //     FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
  //     fe_face_values.reinit(cell, faceID);
  //     this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, c_dof, R, flux); // minus flux for inbound // plus fulx for outbound
  //   }
  // }

}

