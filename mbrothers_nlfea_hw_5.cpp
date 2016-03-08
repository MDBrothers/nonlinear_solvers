#include <iostream>
#include <vector>
#include <exception>
#include "/home/mbrothers/Projects/armadillo/include/armadillo"
//#include <armadillo>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>

#define MAX_ITRS 14

/* File output */
template <typename scalarT>
void writeCSVContents(const std::string filename, std::vector<scalarT>& container, unsigned NCOLS)
{
  try{

      std::ofstream outfile;
    outfile.open( filename, std::ios::app );

    const unsigned ContainerLength= container.size();
    if(not ((ContainerLength % NCOLS)==0)) throw 121;

    if(outfile){
      for(unsigned containerIndex(0); containerIndex<ContainerLength; ){
        for(unsigned lineIndex(0); lineIndex < (NCOLS-1); ++lineIndex){
          outfile << container[containerIndex++] << ", ";
        }
        outfile << container[containerIndex++] << "\n";
      }
    }
    outfile.close();
  }
  catch(int e){
    std::cout << "An exception occured. Nr. " << e << '\n';
  }
}
//template void writeCSVContents<double>( const std::string filename, std::vector<double>& container, unsigned NCOLS);
//template void writeCSVContents<std::string>( const std::string filename, std::vector<std::string>& container, unsigned NCOLS);


/* Per the problem statement, compute the elements of N */
void eval_N(arma::vec::fixed<2> &myN,
                const arma::vec::fixed<2> &myd,
                const double myx)
{
  myN[0] = myx*myd[0]*myd[0]*myd[0] - 2.0*myd[1]*myd[1] + 5.6*myd[0];
  myN[1] = myd[1] - myd[0];
}

/* Compute the elements of the tangent stiffness matrix using centered finite difference. */
void compute_jacobian(arma::mat::fixed<2,2> &myTS,
                               arma::vec::fixed<2> &myd, const double myx,
                               const double myeps)
{
  arma::vec::fixed<2> neg_pert_N, pos_pert_N;
  arma::vec::fixed<2> unperturbed_d;

  unperturbed_d = myd;

  myd[0] += myeps;
  eval_N(pos_pert_N, myd, myx);
  myd[0] = unperturbed_d[0];

  myd[0] -= myeps;
  eval_N(neg_pert_N, myd, myx);
  myd[0] = unperturbed_d[0];

  /* dN1/dd1 and dN2/dd1 */
  myTS.col(0) = (pos_pert_N - neg_pert_N)/(2.0*myeps);

  myd[1] += myeps;
  eval_N(pos_pert_N, myd, myx);
  myd[1] = unperturbed_d[1];

  myd[1] -= myeps;
  eval_N(neg_pert_N, myd, myx);
  myd[1] = unperturbed_d[1];

  /* dN1/dd2 and dN2/dd2 */
  myTS.col(1) = (pos_pert_N - neg_pert_N)/(2.0*myeps);
}

/* eval the res according to the criterion as in the problem statement. */
bool converged(const arma::vec::fixed<2> &myres,
                              const double myeps,
                              const double frst_res_nrm,
                              double &myres_nrm)
{
  myres_nrm = arma::norm(myres, 2);
  return myres_nrm <= (myeps*frst_res_nrm);
}

void save_soln(std::vector< arma::vec::fixed<2> > &myd_vectors,
                  const arma::vec::fixed<2> &myd)
{
  myd_vectors.push_back(myd);
}

void comp_and_store_v_BFGS(std::vector<arma::vec::fixed<2> > &myv_BFGS,
                           const arma::vec::fixed<2> &mydelta_d,
                           const arma::vec::fixed<2> &mydelta_R)
{
  myv_BFGS.push_back(mydelta_d/(arma::dot(mydelta_d, mydelta_R)));
}

void comp_and_store_w_BFGS(std::vector<arma::vec::fixed<2> > &myw_BFGS,
                           const arma::vec::fixed<2> &mydelta_R,
                           const arma::vec::fixed<2> &myR,
                           const double myalpha)
{
  myw_BFGS.push_back(myalpha*myR - mydelta_R);
}

void comp_and_store_alpha_BFGS(std::vector<double> &myalpha_BFGS,
                               const arma::vec::fixed<2> &mydelta_d,
                               const arma::vec::fixed<2> &mydelta_R,
                               const arma::vec::fixed<2> &myR,
                               const double mys)
{
  myalpha_BFGS.push_back(std::sqrt(-mys*arma::dot(mydelta_R, mydelta_d)/
                                        arma::dot(myR, mydelta_d)));
}

int main(){
  /* Basic data structures used for each problem part */
  arma::mat::fixed<2,2> jacobian;
  arma::vec::fixed<2> int_frc, d, delta_d, ext_frc;

  /* For each part of the problem, 'numerical' d values.*/
  std::vector< arma::vec::fixed<2> > d_pure_NR;
  std::vector< arma::vec::fixed<2> > d_moded_NR;
  std::vector< arma::vec::fixed<2> > d_moded_NR_wls;
  std::vector< arma::vec::fixed<2> > d_moded_NR_BFGS;
  std::vector< arma::vec::fixed<2> > d_moded_NR_BFGS_wls;

  /* Num iterations per load step */
  std::vector<int> itrs_pure_NR;
  std::vector<int> itrs_moded_NR;
  std::vector<int> itrs_moded_NR_wls;
  std::vector<int> itrs_moded_NR_BFGS;
  std::vector<int> itrs_moded_NR_BFGS_wls;

  /* The BFGS vectors for part four */
  std::vector< arma::vec::fixed<2> > BFGS_v;
  std::vector< arma::vec::fixed<2> > BFGS_w;
  std::vector< double > BFGS_alpha;

  /* Initialize the reusable data structures to zero and give values to the control constants*/
  jacobian.zeros(); int_frc.zeros(); d.zeros(); delta_d.zeros();
  ext_frc.zeros();
  const int num_load_steps(40);
  const double load_increment(0.25), x(0.19), eps_conv(1.E-4), eps_probe(1.E-8);

  pure:
  /********************************/
  /*Solve with pure Newton Raphson*/
  /********************************/
  std::cout << "\nNew Method: pure Newton Raphson" << std::endl;

  /* Apply each of the load steps and attempt to find d to the specified tolerance. */
  for(int step(0); step < num_load_steps; step++){

    /* Update the load vector */
    ext_frc[0] = (1+step)*load_increment;

    /*  Compute the frst res of the load step*/
    eval_N(int_frc, d, x);
    double frst_res_nrm = arma::norm(ext_frc - int_frc, 2);
    double res_nrm(0.0);
    std::cout << "\nStep: " << step + 1 << " Conv criterion: "
              << frst_res_nrm*eps_conv << std::endl;

    /* Attempt to iteratively solve for d, if we fail, . */
    for(int itr(0); itr < MAX_ITRS; itr++){

      /* Compute the consistent tangent stiffness matrix */
      compute_jacobian(jacobian, d, x, eps_probe);

      /* Solve for the update vector, then update d */
      d += arma::solve(jacobian, ext_frc - int_frc);

      /* Have we reached an acceptable soln? */
      eval_N(int_frc, d, x);
      if(converged(ext_frc - int_frc, eps_conv, frst_res_nrm, res_nrm)){
        save_soln(d_pure_NR, d);
        itrs_pure_NR.push_back(itr + 1);
        std::cout <<"\tconverged: " << itr << " res: " << res_nrm << std::endl;
        break; //Exit inner loop and re-enter the load step loop because we have succeeded.
      }
      std::cout <<"\titr: " << itr << " res: " << res_nrm << std::endl;
      if(res_nrm!=res_nrm)
      { std::cout << "**** NAN in resdiual eval, skipping" << std::endl;
      goto moded;}

      /* The method as applied has broken down, warn the user  */
      if(itr == MAX_ITRS - 1)
      {std::cout << "**** max num ITRS exceeded, skipping" << std::endl;
      goto moded;}
    }
  }

  moded:

  /************************************/
  /*Solve with moded Newton Raphson   */
  /************************************/
  std::cout << "\nNew Method: moded Newton Raphson" << std::endl;
  jacobian.zeros(), int_frc.zeros(), d.zeros(), delta_d.zeros(), ext_frc.zeros();

  /* Apply each of the load steps and attempt to find d to the specified tolerance. */
  for(int step(0); step < num_load_steps; step++){

    /* Update the load vector */
    ext_frc[0] = (1+step)*load_increment;

    /* Compute the frst res of the load step */
    eval_N(int_frc, d, x);
    double frst_res_nrm = arma::norm(ext_frc - int_frc, 2);
    double res_nrm(0.0);
    std::cout << "\nStep: " << step + 1 << " Conv criterion: "
                            << frst_res_nrm*eps_conv << std::endl;

    /* Attempt to iteratively solve for d, if we fail, . */
    for(int itr(0); itr < MAX_ITRS; itr++){

      /* Compute the consistent tangent stiffness matrix but only on the frst itr*/
      if(itr == 0) compute_jacobian(jacobian, d, x, eps_probe);

      /*  Solve for the update vector, then update d */
      d += arma::solve(jacobian, ext_frc - int_frc);

      /* Have we reached an acceptable soln? */
      eval_N(int_frc, d, x);
      if(converged(ext_frc - int_frc, eps_conv, frst_res_nrm, res_nrm)){
        save_soln(d_moded_NR, d);
        itrs_moded_NR.push_back(itr + 1);
        std::cout <<"\tconverged: " << itr << " res: " << res_nrm << std::endl;
        break; //Exit inner loop and re-enter the load step loop because we have succeeded.
      }
      std::cout <<"\titr: " << itr << " res: " << res_nrm << std::endl;
      if(res_nrm!=res_nrm)
      {std::cout << "**** NAN in resdiual eval, skipping" << std::endl;
      goto moded_wls;}

      /* The method as applied has broken down, warn the user  */
      if(itr == MAX_ITRS - 1)
      {std::cout << "**** max num ITRS exceeded, skipping" << std::endl;
      goto moded_wls;}
    }
  }

  moded_wls:

  /*****************************************************/
  /*Solve with moded Newton Raphson with line search*/
  /*****************************************************/
  std::cout << "\nNew Method: moded Newton Raphson with line search" << std::endl;
  jacobian.zeros(), int_frc.zeros(), d.zeros(), delta_d.zeros(), ext_frc.zeros();

  /* Apply each of the load steps and attempt to find d to the specified tolerance. */
  for(int step(0); step < num_load_steps; step++){

    /* Update the load vector */
    ext_frc[0] = (1+step)*load_increment;

    /* Compute the frst res of the load step */
    eval_N(int_frc, d, x);
    double frst_res_nrm = arma::norm(ext_frc - int_frc, 2);
    double res_nrm(0.0);
    std::cout << "\nStep: " << step + 1 << " Conv criterion: "
              << frst_res_nrm*eps_conv << std::endl;

    /* Attempt to iteratively solve for d, if we fail, . */
    for(int itr(0); itr < MAX_ITRS; itr++){

      /* Compute the consistent tangent stiffness matrix but only on the frst itertation */
      if(itr == 0) compute_jacobian(jacobian, d, x, eps_probe);

      /* Solve for the update vector, then initialize the line search */
      delta_d = arma::solve(jacobian, ext_frc - int_frc);
      double srch_prm(1.0), G_zero(0);
      bool ln_srch_win(false);
      G_zero = arma::dot(delta_d, ext_frc - int_frc);

      /* Perform the line search, exit if failed */
      for(int search_it(0); search_it < 5; search_it++){
        eval_N(int_frc, d + srch_prm*delta_d, x);
        if(std::abs(arma::dot(delta_d, ext_frc - int_frc)) <
           std::abs(.5*G_zero)){
          d += delta_d*srch_prm; ln_srch_win = true; break;
        }
        else srch_prm*=1.0/std::sqrt(2.0); //Always try a smaller step.
      }
      if(ln_srch_win == false)
      {std::cout << "**** line search failed, skipping" << std::endl;
      goto moded_bfgs;}

      /* Have we reached an acceptable soln? */
      eval_N(int_frc, d, x);
      if(converged(ext_frc - int_frc, eps_conv, frst_res_nrm, res_nrm)){
        save_soln(d_moded_NR_wls, d);
        itrs_moded_NR_wls.push_back(itr + 1);
        std::cout <<"\tconverged: " << itr << " res: " << res_nrm << std::endl;
        break; //Exit inner loop and re-enter the load step loop because we have succeeded.
      }
      std::cout <<"\titr: " << itr << " res: " << res_nrm << std::endl;
      if(res_nrm!=res_nrm)
      {std::cout << "**** NAN in resdiual eval, skipping" << std::endl;
      goto moded_bfgs;}

      /* The method as applied has broken down, warn the user  */
      if(itr == MAX_ITRS - 1)
      {std::cout << "**** max num ITRS exceeded, skipping" << std::endl;
      goto moded_bfgs;}
    }
  }

  moded_bfgs:

  /**********************************************/
  /*Solve with moded Newton Raphson with BFGS*/
  /**********************************************/
  jacobian.zeros(), int_frc.zeros(), d.zeros(), delta_d.zeros(), ext_frc.zeros();
  std::cout << "\nNew Method: moded Newton Raphson with BFGS" << std::endl;

  /* Some variables needed for BFGS */
  arma::vec::fixed<2> delta_R;
  arma::vec::fixed<2> previous_R;
  arma::vec::fixed<2> current_R;
  arma::mat::fixed<2,2> lh_BFGS_matrix, rh_BFGS_matrix;

  /* Apply each of the load steps and attempt to find d to the specified tolerance. */
  for(int step(0); step < num_load_steps; step++){

    /* Clear previous accumulated BFGS update matricesy */
    delta_R.zeros();
    previous_R.zeros();
    current_R.zeros();
    rh_BFGS_matrix = arma::eye< arma::mat >(2,2);
    lh_BFGS_matrix = arma::eye< arma::mat >(2,2);
    BFGS_v.clear();
    BFGS_w.clear();
    BFGS_alpha.clear();

    /* Update the load vector */
    ext_frc[0] = (1+step)*load_increment;

    /* Compute the frst res of the load step */
    eval_N(int_frc, d, x);
    double frst_res_nrm = arma::norm(ext_frc - int_frc, 2);
    double res_nrm(0.0);
    std::cout << "\nStep: " << step + 1 << " Conv criterion: " << frst_res_nrm*eps_conv << std::endl;

    /* Attempt to iteratively solve for d, if we fail, . */
    for(int itr(0); itr < MAX_ITRS; itr++){

      /* Compute the res for no update in d */
      eval_N(int_frc, d, x);
      previous_R = ext_frc - int_frc;

      /* Compute the consistent tangent stiffness matrix but only on the frst itr*/
      if(itr == 0){
        compute_jacobian(jacobian, d, x, eps_probe);

        /* Invert the jacobian and solve for delta_d */
        jacobian = arma::inv(jacobian);
        delta_d = jacobian*(ext_frc - int_frc);
        d += delta_d;
      }
      else{

        /* solve for delta_d, then update d */
        delta_d = lh_BFGS_matrix*jacobian*rh_BFGS_matrix*current_R;
        d += delta_d;
      }

      /* Store the update information */
      eval_N(int_frc, d, x);
      current_R = ext_frc - int_frc;
      delta_R = current_R - previous_R;

      /* compute the BFGS update vectors */
      comp_and_store_alpha_BFGS(BFGS_alpha, delta_d, delta_R, previous_R, 1.0);
      comp_and_store_v_BFGS(BFGS_v, delta_d, delta_R);
      comp_and_store_w_BFGS(BFGS_w, delta_R, previous_R, BFGS_alpha.back());

      /* Update the BFGS update matrices */
      lh_BFGS_matrix = (arma::eye< arma::mat >(2,2) +
                        BFGS_v.back()*(BFGS_w.back().t()) )*lh_BFGS_matrix;
      rh_BFGS_matrix = rh_BFGS_matrix*(arma::eye< arma::mat >(2,2) +
                       BFGS_w.back()*(BFGS_v.back().t()) );

      /* Test for Conv */
      if(converged(ext_frc - int_frc, eps_conv, frst_res_nrm, res_nrm)){
        save_soln(d_moded_NR_BFGS, d);
        itrs_moded_NR_BFGS.push_back(itr + 1);
        std::cout <<"\tconverged: " << itr << " res: " << res_nrm << std::endl;
        break; //Exit inner loop and re-enter the load step loop because we have succeeded.
      }
      std::cout <<"\titr: " << itr << " res: " << res_nrm << std::endl;
      if(res_nrm!=res_nrm)
      {std::cout << "**** NAN in resdiual eval, skipping" << std::endl;
      goto moded_bfgs_wls;}

      /* The method as applied has broken down, warn the user  */
      if(itr == MAX_ITRS - 1)
      {std::cout << "**** max num ITRS exceeded, skipping" << std::endl;
      goto moded_bfgs_wls;}
    }
  }

  moded_bfgs_wls:

  /**************************************************************/
  /*Solve with moded Newton Raphson with BFGS and line search*/
  /**************************************************************/
  jacobian.zeros(), int_frc.zeros(), d.zeros(), delta_d.zeros(), ext_frc.zeros();
  std::cout << "\nNew Method: moded Newton Raphson with BFGS and line search" << std::endl;

  /* Apply each of the load steps and attempt to find d to the specified tolerance. */
  for(int step(0); step < num_load_steps; step++){

    /* Clear previous accumulated BFGS update matricesy */
    delta_R.zeros();
    previous_R.zeros();
    current_R.zeros();
    rh_BFGS_matrix = arma::eye< arma::mat >(2,2);
    lh_BFGS_matrix = arma::eye< arma::mat >(2,2);
    BFGS_v.clear();
    BFGS_w.clear();
    BFGS_alpha.clear();

    /* Update the load vector */
    ext_frc[0] = (1+step)*load_increment;

    /* Compute the frst res of the load step */
    eval_N(int_frc, d, x);
    double frst_res_nrm = arma::norm(ext_frc - int_frc, 2);
    double res_nrm(0.0);
    std::cout << "\nStep: " << step + 1 << " Conv criterion: "
              << frst_res_nrm*eps_conv << std::endl;

    /* Attempt to iteratively solve for d,*/
    for(int itr(0); itr < MAX_ITRS; itr++){

      /* Compute the res for no update in d */
      eval_N(int_frc, d, x);
      previous_R = ext_frc - int_frc;

      /* Compute the consistent tangent stiffness matrix but only on the frst itr*/
      if(itr == 0){
        compute_jacobian(jacobian, d, x, eps_probe);

        /* Invert the jacobian and solve for delta_d */
        jacobian = arma::inv(jacobian);
        delta_d = jacobian*(ext_frc - int_frc);
      }
      else{
        /* solve for delta_d, then update d */
        delta_d = lh_BFGS_matrix*jacobian*rh_BFGS_matrix*current_R;
      }

      /* Initialize the line search */
      double srch_prm(1.0), G_zero(0);
      bool ln_srch_win(false);
      G_zero = arma::dot(delta_d, previous_R);

      /* Perform the line search, exit if failed */
      for(int search_it(0); search_it < 5; search_it++){
        eval_N(int_frc, d + srch_prm*delta_d, x);
        if(std::abs(arma::dot(delta_d, ext_frc - int_frc)) <
           std::abs(.5*G_zero)){
          d += delta_d*srch_prm; ln_srch_win = true; break;
        }
        else srch_prm*=1.0/std::sqrt(2.0); //Always try a smaller step.
      }
      if(ln_srch_win == false)
      {std::cout << "**** line search failed, skipping" << std::endl; goto end;}

      /* Store the update information */
      current_R = ext_frc - int_frc;
      delta_R = current_R - previous_R;

      /* compute the BFGS update vectors */
      comp_and_store_alpha_BFGS(BFGS_alpha, delta_d, delta_R, previous_R, srch_prm);
      comp_and_store_v_BFGS(BFGS_v, delta_d, delta_R);
      comp_and_store_w_BFGS(BFGS_w, delta_R, previous_R, BFGS_alpha.back());

      /* Update the BFGS update matrices */
      lh_BFGS_matrix = (arma::eye< arma::mat >(2,2) +
                        BFGS_v.back()*(BFGS_w.back().t()) )*lh_BFGS_matrix;
      rh_BFGS_matrix = rh_BFGS_matrix*(arma::eye< arma::mat >(2,2) +
                       BFGS_w.back()*(BFGS_v.back().t()) );

      /* Test for Conv */
      if(converged(ext_frc - int_frc, eps_conv, frst_res_nrm, res_nrm)){
        save_soln(d_moded_NR_BFGS_wls, d);
        itrs_moded_NR_BFGS_wls.push_back(itr + 1);
        std::cout <<"\tconverged: " << itr << " res: " << res_nrm << std::endl;
        break; //Exit inner loop and re-enter the load step loop because we have succeeded.
      }
      std::cout <<"\titr: " << itr << " res: " << res_nrm << std::endl;
      if(res_nrm!=res_nrm)
      {std::cout << "**** NAN in resdiual eval, skipping" << std::endl;
      goto end;}

      /* The method as applied has broken down, warn the user  */
      if(itr == MAX_ITRS - 1)
      {std::cout << "**** max num ITRS exceeded, skipping" << std::endl;
      goto end;}
    }
  }

  end:

  /*output the converged solutions*/

  std::vector<std::string> labels;
  labels.push_back("d_NR");
  labels.push_back("N_NR");
  writeCSVContents<std::string>("pure_nr_x19.txt", labels, 2);
  for(const auto & soln : d_pure_NR){
    std::vector<double> record;
    record.push_back( soln[0]);
    d[0] = soln[0]; d[1] = soln[0];
    eval_N(int_frc, d, x);
    record.push_back(int_frc[0]);
    writeCSVContents<double>("pure_nr_x19.txt",record,2);
  }

  labels.clear();
  labels.push_back("num iters pure nr");
  writeCSVContents<std::string>("pure_nr_x19_conv.txt", labels, 1);
  writeCSVContents<int>("pure_nr_x19_conv.txt", itrs_pure_NR, 1);

  labels.clear();
  labels.push_back("d_moded_NR");
  labels.push_back("N_moded_NR");
  writeCSVContents<std::string>("moded_nr_x19.txt", labels, 2);
  for(const auto & soln : d_moded_NR){
    std::vector<double> record;
    record.push_back( soln[0]);
    d[0] = soln[0]; d[1] = soln[0];
    eval_N(int_frc, d, x);
    record.push_back(int_frc[0]);
    writeCSVContents<double>("moded_nr_x19.txt",record,2);
  }

  labels.clear();
  labels.push_back("num iters moded nr");
  writeCSVContents<std::string>("moded_nr_x19_conv.txt", labels, 1);
  writeCSVContents<int>("moded_nr_x19_conv.txt", itrs_moded_NR, 1);

  labels.clear();
  labels.push_back("d_moded_NR_wls");
  labels.push_back("N_moded_NR_wls");
  writeCSVContents<std::string>("moded_nr_wls_x19.txt", labels, 2);
  for(const auto & soln : d_moded_NR_wls){
    std::vector<double> record;
    record.push_back( soln[0]);
    d[0] = soln[0]; d[1] = soln[0];
    eval_N(int_frc, d, x);
    record.push_back(int_frc[0]);
    writeCSVContents<double>("moded_nr_wls_x19.txt",record,2);
  }

  labels.clear();
  labels.push_back("num iters moded nr wls");
  writeCSVContents<std::string>("moded_nr_wls_x19_conv.txt", labels, 1);
  writeCSVContents<int>("moded_nr_wls_x19_conv.txt", itrs_moded_NR_wls, 1);

  labels.clear();
  labels.push_back("d_moded_NR_BFGS");
  labels.push_back("N_moded_NR_BFGS");
  writeCSVContents<std::string>("moded_nr_bfgs_x19.txt", labels, 2);
  for(const auto & soln : d_moded_NR_BFGS){
    std::vector<double> record;
    record.push_back( soln[0]);
    d[0] = soln[0]; d[1] = soln[0];
    eval_N(int_frc, d, x);
    record.push_back(int_frc[0]);
    writeCSVContents<double>("moded_nr_bfgs_x19.txt",record,2);
  }

  labels.clear();
  labels.push_back("num iters moded nr bfgs");
  writeCSVContents<std::string>("moded_nr_bfgs_x19_conv.txt", labels, 1);
  writeCSVContents<int>("moded_nr_bfgs_x19_conv.txt", itrs_moded_NR_BFGS, 1);

  labels.clear();
  labels.push_back("d_moded_NR_BFGS_wls");
  labels.push_back("N_moded_NR_BFGS_wls");
  writeCSVContents<std::string>("moded_nr_bfgs_wls_x19.txt", labels, 2);
  for(const auto & soln : d_moded_NR_BFGS_wls){
    std::vector<double> record;
    record.push_back( soln[0]);
    d[0] = soln[0]; d[1] = soln[0];
    eval_N(int_frc, d, x);
    record.push_back(int_frc[0]);
    writeCSVContents<double>("moded_nr_bfgs_wls_x19.txt",record,2);
  }

  labels.clear();
  labels.push_back("num iters moded nr bfgs wls");
  writeCSVContents<std::string>("moded_nr_bfgs_wls_x19_conv.txt", labels, 1);
  writeCSVContents<int>("moded_nr_bfgs_wls_x19_conv.txt", itrs_moded_NR_BFGS_wls, 1);

  return 0;
}
