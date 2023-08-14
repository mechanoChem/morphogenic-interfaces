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
#include <random> // For random number generation

template <int dim>
class CahnHilliard : public mechanoChemFEM<dim>
{
public:
  CahnHilliard();
  // this is a overloaded function
  void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
                    const FEValues<dim> &fe_values,
                    Table<1, Sacado::Fad::DFad<double>> &R,
                    Table<1, Sacado::Fad::DFad<double>> &ULocal,
                    Table<1, double> &ULocalConv);

  void apply_initial_condition();
  // void impose_constraints();
  void solve_ibvp();
  void output_results();

  ParameterHandler *params;
  ConstraintMatrix *constraints;

  int c1_dof = 0, mu1_dof = 1, c2_dof = 2, mu2_dof = 3,  out1_dof = 4, out2_dof = 5;
  
  int iter_count = 0;
  bool increase_ts = false;

  // std::vector<double> x, y, r, v;

  Vector<double> theta_1;
  Vector<double> theta_2;
  Vector<double> theta_3;
  Vector<double> theta_4;
};

template <int dim>
CahnHilliard<dim>::CahnHilliard()
{
  // This let you use one params to get all parameters pre-defined in the mechanoChemFEM
  params = this->params_mechanoChemFEM;
  params->enter_subsection("Concentration");
  //
  params->declare_entry("xw", "1.0", Patterns::FileName());
  params->declare_entry("yw", "1.0", Patterns::FileName());
  params->declare_entry("rw", "1.0", Patterns::FileName());
  params->declare_entry("vw", "1.0", Patterns::FileName());
  //
  params->declare_entry("c1_ini", "0.25", Patterns::Double());
  params->declare_entry("c2_ini", "0.25", Patterns::Double());
  params->declare_entry("c1_per", "0.01", Patterns::Double());
  params->declare_entry("c2_per", "0.01", Patterns::Double());
  //
  params->declare_entry("mobility_1","0.1",Patterns::Double() );
  params->declare_entry("mobility_2","0.1",Patterns::Double() );
  params->declare_entry("kappa_1","1.0",Patterns::Double() );
  params->declare_entry("kappa_2","1.0",Patterns::Double() );
   // 
  params->declare_entry("flux_time","100",Patterns::Double() );
  params->declare_entry("flux_amp","1.0",Patterns::Double() );
  params->leave_subsection();

  // Declare the parameters before load it
  this->load_parameters("parameters.prm");

  // define main fields from parameter file.
  this->define_primary_fields();
  // Set up the ibvp.
  this->init_ibvp();

  // -----------------------------------------------------------------------
  // add constraints and BC
  constraints = this->constraints_mechanoChemFEM;
  constraints->clear();
  DoFTools::make_hanging_node_constraints(this->dof_handler, *constraints);
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
          if (ck == out1_dof || ck == out2_dof)
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
  this->pcout << "===========  imposing constraints ===========" << std::endl;
  //-----------------------------------------------------------------------

  theta_1.reinit(this->triangulation.n_active_cells());
  theta_2.reinit(this->triangulation.n_active_cells());
  theta_3.reinit(this->triangulation.n_active_cells());
  theta_4.reinit(this->triangulation.n_active_cells());
}

template <int dim>
void CahnHilliard<dim>::apply_initial_condition()
{
  this->pcout << "=========== applying initial condition (new) ===========" << std::endl;
  // params->enter_subsection("Geometry");
  // double x_min = params->get_double("x_min");
  // double x_max = params->get_double("x_max");
  // double y_min = params->get_double("y_min");
  // double y_max = params->get_double("y_max");
  // params->leave_subsection();
  //
  params->enter_subsection("Concentration");
  double c1_ini = params->get_double("c1_ini");
  double c2_ini = params->get_double("c2_ini");
  double c1_per = params->get_double("c1_per");
  double c2_per = params->get_double("c2_per");
  params->leave_subsection();

  unsigned int seed = 12345;
  std::mt19937 generator(seed);
  std::uniform_real_distribution<double> distribution(-0.5, 0.5);

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
          
          this->solution_prev(local_dof_indices[i]) = c1_ini + c1_per * distribution(generator);
        if (ck == mu1_dof)
          this->solution_prev(local_dof_indices[i]) = 0.0;
        if (ck == c2_dof)
          this->solution_prev(local_dof_indices[i]) = c2_ini + c2_per * distribution(generator);
        if (ck == mu2_dof)
          this->solution_prev(local_dof_indices[i]) = 0.0;
        if (ck == out1_dof)
          this->solution_prev(local_dof_indices[i]) = 0.0;
        if (ck == out2_dof)
          this->solution_prev(local_dof_indices[i]) = 0.0;

      } // dofs_per_cell
    }   // this_mpi_process
  }     // cell

  this->solution_prev.compress(VectorOperation::insert); // potentially for parallel computing purpose.
  this->solution = this->solution_prev;

  
}

// template <int dim>
// void CahnHilliard<dim>::impose_constraints()
// {
//   constraints->clear();
//   DoFTools::make_hanging_node_constraints(this->dof_handler, *constraints);
//   //-------------------
//   // add constraints and BC
//   //-------------------
//   {
//     hp::FEValues<dim> hp_fe_values(this->fe_collection, this->q_collection, update_values | update_quadrature_points);
//     typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
//     for (; cell != endc; ++cell)
//     {
//       if (cell->subdomain_id() == this->this_mpi_process)
//       { // parallel computing.
//         // int cell_id = cell->active_cell_index();
//         hp_fe_values.reinit(cell);
//         const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

//         const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
//         std::vector<unsigned int> local_dof_indices(dofs_per_cell);
//         cell->get_dof_indices(local_dof_indices);
//         for (unsigned int i = 0; i < dofs_per_cell; ++i)
//         {
//           const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first;
//           // std::cout << cell_id << " ck " << ck << std::endl;
//           if (ck == out1_dof || ck == out2_dof)
//           {
//             auto globalDOF = local_dof_indices[i];
//             constraints->add_line(globalDOF);
//             constraints->set_inhomogeneity(globalDOF, 0.0);
//           }
//         }
//       }
//     }
//   }
//   constraints->close();
//   std::cout << "===========  end of constraint =========== " << std::endl;
// }

template <int dim>
void CahnHilliard<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim> &fe_values, Table<1, Sacado::Fad::DFad<double>> &R, Table<1, Sacado::Fad::DFad<double>> &ULocal, Table<1, double> &ULocalConv)
{

  params->enter_subsection("Geometry");
  // double x_min = params->get_double("x_min");
  // double x_max = params->get_double("x_max");
  // double y_min = params->get_double("y_min");
  // double y_max = params->get_double("y_max");
  double z_min = params->get_double("z_min");
  // double z_max = params->get_double("z_max");
  // int num_elem_x = params->get_integer("num_elem_x");
  // int num_elem_y = params->get_integer("num_elem_y");
  params->leave_subsection();

  // evaluate primary fields
  params->enter_subsection("Concentration");
  //
  std::vector<double> x = get_double_vector(params->get("xw"));
  std::vector<double> y = get_double_vector(params->get("yw"));
  std::vector<double> r = get_double_vector(params->get("rw"));
  std::vector<double> v = get_double_vector(params->get("vw"));
  //
  double M1 = params->get_double("mobility_1");
	double M2 = params->get_double("mobility_2");
	double k1 = params->get_double("kappa_1");
	double k2 = params->get_double("kappa_2");
  // double test = params->get_double("kappa_2");
  //
  double flux_time = params->get_double("flux_time");
  double flux_amp = params->get_double("flux_amp");
  params->leave_subsection();

  // std::vector<double> x = {0.1, 0.90, 0.10, 0.45};
  // std::vector<double> y = {0.10, 0.10, 0.90, 0.45};
  // std::vector<double> v = {-1.0, -1.0, -1.0, 0.2};
  // std::vector<double> r = {8.0, 8.0, 8.0, 1.0};
  // double a = 1.0;
  // double m = -1.0; //= min v_val

  double T = 300; // K

  /**************************************
  * Define chemistry
  *************************************/
  double a0 = -1.51629e+10 / 1e27 * 6.242e18; // ev/nm3
  double b0 = 1.95270e+10 / 1e27 * 6.242e18;  // ev/nm3
  double c0 = -4.36414e+09 / 1e27 * 6.242e18; // ev/nm3

  // double k = 1.380649e-23 * 6.242e18;        // m2 kg s-2 K-1 = J/K => eV/K;
  double k = 8.617333262e-5;     // eV/K
  double V = 0.0201207;          // nm^3   = 20.1207e-30 * 1e27;

  double d = 100;
  // --------------------------------------------------------------------

 	unsigned int n_q_points= fe_values.n_quadrature_points;

  //define fields
	dealii::Table<1,double>  c1_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > c1(n_q_points), mu1(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > c1_grad(n_q_points, dim), mu1_grad(n_q_points, dim);
	
	dealii::Table<1,double>  c2_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > c2(n_q_points), mu2(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > c2_grad(n_q_points, dim), mu2_grad(n_q_points, dim);
  //evaluate fields
	evaluateScalarFunction<double,dim>(fe_values, c1_dof, ULocalConv, c1_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c1_dof, ULocal, c1);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c1_dof, ULocal, c1_grad);
  //
	evaluateScalarFunction<double,dim>(fe_values, c2_dof, ULocalConv, c2_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c2_dof, ULocal, c2);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c2_dof, ULocal, c2_grad);
	//
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, mu1_dof, ULocal, mu1);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, mu1_dof, ULocal, mu1_grad);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, mu2_dof, ULocal, mu2);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, mu2_dof, ULocal, mu2_grad);
	
	//evaluate diffusion and reaction term
	dealii::Table<1,Sacado::Fad::DFad<double> > rhs_mu1(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_c1(n_q_points, dim), kappa_c1_grad(n_q_points, dim);
	dealii::Table<1,Sacado::Fad::DFad<double> > rhs_mu2(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_c2(n_q_points, dim), kappa_c2_grad(n_q_points, dim);
  Sacado::Fad::DFad<double> ff_1, ff_2, f_c1, f_c2, F_c1, F_c2;

  // double d = .4;
  // double s = .7;
	for(unsigned int q=0; q < n_q_points; q++){

    copy_to_q(j_c1, -M1 * copy_from_q(mu1_grad, q), q);  //j_c1[q][i] = -M_1 * c_1_grad[q][i], for i=1:dim
    copy_to_q(j_c2, -M2 * copy_from_q(mu2_grad, q), q);
    copy_to_q(kappa_c1_grad, k1 * copy_from_q(c1_grad, q), q);
    copy_to_q(kappa_c2_grad, k2 * copy_from_q(c2_grad, q), q);
		
    // F_c1 = 2 * (c1[q] - .5) * std::pow(c2[q] - .5, 2);
    // F_c2 = 2 * (c2[q] - .5) * std::pow(c1[q] - .5, 2);
    // for (unsigned int i = 0; i < x.size(); i++){
    //   F_c1 += -r[i] * v[i] * exp(-r[i]*(std::pow(c1[q] - x[i],2) + std::pow(c2[q] - y[i],2))) * (2*c1[q] - 2*x[i]);
    //   F_c2 += -r[i] * v[i] * exp(-r[i]*(std::pow(c1[q] - x[i],2) + std::pow(c2[q] - y[i],2))) * (2*c2[q] - 2*y[i]);
    // }
    // F_c1 = a * F_c1;
		// F_c2 = a * F_c2;

    // F_c1 = 6*d/std::pow(s,4)*(c1[q]*c1[q]+c2[q]*c2[q])*c1[q]-6*d/std::pow(s,3)*c1[q]*c2[q]-3*d/std::pow(s,2)*c1[q];
    // F_c2 = 6*d/std::pow(s,4)*(c1[q]*c1[q]+c2[q]*c2[q])*c2[q]+3*d/std::pow(s,3)*c2[q]*c2[q]-3*d/std::pow(s,2)*c2[q];
		
    f_c1 = 0.865 * (a0 + b0 + 2 * c0 * c1[q] + (k * T/ V * d) * (exp(d * (c1[q] - 1)) - exp(-d * c1[q])));   // eV/nm^3 
    f_c1 += -(20.0*c2[q]*(c1[q] + 2*c2[q] - 1.0)) / std::pow(c1[q] - 1.0, 3); // eV/nm^3 
   
    // f_c2 = 40.0 * c2[q] - 20.0; // eV/nm^3  
    f_c2 = (20.0*(c1[q] + 2*c2[q] - 1.0))/ std::pow(c1[q] - 1.0, 2); // eV/nm^3 
    
    F_c1 = f_c1;
    F_c2 = f_c2; 


    rhs_mu1[q] = F_c1 - mu1[q];
    rhs_mu2[q] = F_c2 - mu2[q];
	}
	
	//call residual functions
	this->ResidualEq.residualForDiffusionEq(fe_values, c1_dof, R, c1, c1_conv, j_c1);
  this->ResidualEq.residualForDiffusionEq(fe_values, c2_dof, R, c2, c2_conv, j_c2);

	this->ResidualEq.residualForPoissonEq(fe_values, mu1_dof, R, kappa_c1_grad, rhs_mu1);
	this->ResidualEq.residualForPoissonEq(fe_values, mu2_dof, R, kappa_c2_grad, rhs_mu2);
	
  	//calculate features
  double N00, N10, N01;
	dealii::Table<1,Sacado::Fad::DFad<double> > theta1(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > theta2(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > theta3(n_q_points);
  dealii::Table<1,Sacado::Fad::DFad<double> > theta4(n_q_points);
	// dealii::Table<1,Sacado::Fad::DFad<double> > gamma(n_q_points);
	for(unsigned int q=0; q<n_q_points;q++){
		if      (c1[q] <= .5 and c2[q] <= .5) theta1[q]=1;
		else if (c1[q] >  .5 and c2[q] <= .5) theta2[q]=1;
		else if (c1[q] <= .5 and c2[q] >  .5) theta3[q]=1;

    N10 = (c1[q].val() - x[0]) / (x[1] - x[0]);
    N01 = (c2[q].val() - y[0]) / (y[2] - y[0]);
    N00 = 1.0 - N10 - N01;
        
    theta4[q] = N10 * (-1.) +  N00 * (0.) + N01 * (1.);

		// gamma[q]=0.5*(c1_grad[q][0]*c1_grad[q][0]+c1_grad[q][1]*c1_grad[q][1]);
	}
  int cell_id = cell->active_cell_index();
  this->theta_1[cell_id] = this->ResidualEq.volumeIntegration(fe_values, theta1);
  this->theta_2[cell_id] = this->ResidualEq.volumeIntegration(fe_values, theta2);
  this->theta_3[cell_id] = this->ResidualEq.volumeIntegration(fe_values, theta3);
  this->theta_4[cell_id] = this->ResidualEq.volumeIntegration(fe_values, theta4);
	// local_features[0]+=this->ResidualEq.volumeIntegration(fe_values, theta1);
	// local_features[1]+=this->ResidualEq.volumeIntegration(fe_values, theta2);
	// local_features[2]+=this->ResidualEq.volumeIntegration(fe_values, theta3);
	// local_features[3]+=this->ResidualEq.volumeIntegration(fe_values, gamma);


  // -------------------
  // BC
  // -------------------
  double flux = 0.0;
  if (this->current_time > flux_time)
  {
     
    for (unsigned int faceID = 0; faceID < 2 * dim; faceID++)
    {
      if (cell->face(faceID)->at_boundary() == true)
      {
        // this->pcout << c1_avg << c2_avg << std::endl;
         
        if (cell->face(faceID)->center()[2] <= (z_min + 1e-3)) // bottom boundary 
        {
          
          // if (flux_amp > 0)  // adding Li
          // {
            flux = flux_amp;
          // }
          // else if (flux_amp < 0)) // removing Li (only form Li rich area)
          // {
          //   flux = flux_amp;
          // } 
        }
      FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
      fe_face_values.reinit(cell, faceID);
      this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, c1_dof, R, flux); // minus flux for inbound // plus fulx for outbound
      }
    }
  }

}

/**************************************
 * Implement adaptive time stepping
 *************************************/
template <int dim>
void CahnHilliard<dim>::solve_ibvp()
{

  if (increase_ts)
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

  // -------------------------------------------------------

  params->enter_subsection("Concentration");
  std::vector<double> x = get_double_vector(params->get("xw"));
  std::vector<double> y = get_double_vector(params->get("yw"));
  std::vector<double> r = get_double_vector(params->get("rw"));
  std::vector<double> v = get_double_vector(params->get("vw"));
  params->leave_subsection();

  double c1_value, c2_value; 
  // double N00, N10, N01;

  dealii::Vector<double> loc_U(this->solution); // solution at all nodes at current time (on this processor ?)

  const unsigned int vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
  // typename hp::DoFHandler<dim>::active_cell_iterator cell = hpFEM<dim>::dof_handler.begin_active(), endc = hpFEM<dim>::dof_handler.end();
  for (; cell != endc; ++cell)
  {
    if (cell->is_locally_owned()) 
    {
        // const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        // const unsigned int dofs_per_vertex = cell->get_fe().dofs_per_vertex;

      for (unsigned int i = 0; i < vertices_per_cell; ++i)
      {
        c1_value = loc_U(cell->vertex_dof_index(i, c1_dof, cell->active_fe_index()));
        c2_value = loc_U(cell->vertex_dof_index(i, c2_dof, cell->active_fe_index()));
        
        // if      (c1_value <= .5 && c2_value <= .5) theta1[q]=1;
        // else if (c1_value >  .5 && c2_value <= .5) theta2[q]=1;
        // else if (c1_value <= .5 && c2_value >  .5) theta3[q]=1;
      
        // N10 = (c1_value - x[0]) / (x[1] - x[0]);
        // N01 = (c2_value - y[0]) / (y[2] - y[0]);
        // N00 = 1.0 - N10 - N01;
        this->solution(cell->vertex_dof_index(i, out1_dof , cell->active_fe_index())) = 1.0 - c1_value - c2_value; //N10 * (-1.) +  N00 * (0.) + N01 * (1.);
      }
    }
  }
  this->solution.compress(VectorOperation::insert); 
  // -------------------------------------------------------

  this->solution_prev = this->solution;
}

template <int dim>
void CahnHilliard<dim>::output_results()
{
  // Vector<float> subdomain_id(this->triangulation.n_active_cells());
  // typename hp::DoFHandler<dim>::active_cell_iterator elem = this->dof_handler.begin_active(), endc = this->dof_handler.end();
  // unsigned int _j = 0;
  // for (; elem != endc; ++elem)
  // {
  //   subdomain_id(_j++) = elem->subdomain_id();
  // }
  // std::cout << " after elem id " << std::endl;
  Vector<double> _theta_1(this->triangulation.n_active_cells());
  Vector<double> _theta_2(this->triangulation.n_active_cells());
  Vector<double> _theta_3(this->triangulation.n_active_cells());
  Vector<double> _theta_4(this->triangulation.n_active_cells());

  Utilities::MPI::sum(theta_1, MPI_COMM_WORLD, _theta_1);
  Utilities::MPI::sum(theta_2, MPI_COMM_WORLD, _theta_2);
  Utilities::MPI::sum(theta_3, MPI_COMM_WORLD, _theta_3);
  Utilities::MPI::sum(theta_4, MPI_COMM_WORLD, _theta_4);

  // write vtk and snapshot for solution
  if (this->save_output)
  {
    std::string output_path = this->output_directory + "output-" + std::to_string(this->current_increment) + ".vtk";
    this->FEMdata_out.clear_data_vectors();
    this->FEMdata_out.data_out.add_data_vector(_theta_1, "t1");
    this->FEMdata_out.data_out.add_data_vector(_theta_2, "t2");
    this->FEMdata_out.data_out.add_data_vector(_theta_3, "t3");
    this->FEMdata_out.data_out.add_data_vector(_theta_4, "t4");

    if (this->current_increment % this->skip_output == 0)
      this->FEMdata_out.write_vtk(this->solution_prev, output_path);
  }
  if (this->save_snapshot)
  {
    std::string snapshot_path = this->snapshot_directory + "snapshot-" + std::to_string(this->current_increment + this->off_output_index) + ".dat";
    this->pcout << " save to " << snapshot_path << std::endl;
    this->FEMdata_out.create_vector_snapshot(this->solution, snapshot_path);
  }

  // this->pcout
  //     << " total_phi_e " << this->total_phi_e
  //     << " total_phi_b " << this->total_phi_b
  //     << " total_phi_i " << this->total_phi_i
  //     << std::endl;
}