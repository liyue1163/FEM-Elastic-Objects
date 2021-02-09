#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output:
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H) {

   unsigned int itr_cnt = 0;
    Eigen::VectorXd qdot_i = x0;
    g(tmp_g, qdot_i);
    double energy_grad_norm = tmp_g.norm();
    double tol = 1e-2;

    while (energy_grad_norm >= 0.0 && itr_cnt < maxSteps) {
        g(tmp_g, qdot_i);
        H(tmp_H, qdot_i);

        // GiCGSTAB occasionally causes seg fault error.
        // Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(tmp_H);
        Eigen::VectorXd b = -1.0 * tmp_g;

        Eigen::VectorXd d =  solver.solve(b);
        // line search
        double alpha = 1.0;
        double p = 0.5;
        double c = 1e-8;
        double cur_energy = f(qdot_i + alpha * d);
        double threshold = f(qdot_i) + c * (d.transpose() * tmp_g)(0);

        // perform line search
        while (cur_energy > threshold && alpha >= 0.001) {
            alpha = p * alpha;
            cur_energy = f(qdot_i + alpha * d);
        }

        qdot_i = qdot_i + alpha * d;
        g(tmp_g, qdot_i);

        energy_grad_norm = tmp_g.norm();
        itr_cnt = itr_cnt + 1;
    }

    x0 = qdot_i;

    return 0.0;

   // double p = 0.5;
	// double c = std::pow(10.0, -8.0);
	// double tolerence = std::pow(10.0, -1.0);
	// //double tolerence = 0;
	// int step = 0;
	// Eigen::VectorXd x(x0.rows()); //construct x vector to keep track of the current x
	// x = x0;
	// g(tmp_g, x0); //modified
	// while ((step < maxSteps) && (tmp_g.norm() > tolerence)) {
	// 	//calculate hessian and gradient at current point
	// 	g(tmp_g, x);
	// 	H(tmp_H, x);
	// 	//construct a negative g vector
	// 	Eigen::VectorXd neg_tmp_g(tmp_g.rows());
   //
	// 	neg_tmp_g = -1.0 * tmp_g;
	// 	Eigen::VectorXd d; //d is search direction
	// 	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	// 	solver.compute(tmp_H);
	// 	d = solver.solve(neg_tmp_g);
   //
	// 	//choose alpha using line search
	// 	double alpha = 1;//maximum value of alpha is 1
	// 	while (!((f(x + alpha * d) <= f(x) + c * d.transpose() * tmp_g) || (alpha < tolerence))) { //modified
	// 		alpha *= p;
	// 	}
   //
	// 	x += alpha * d;
	// 	step += 1;
   //
	// }
	// x0 = x; //updated x0
	// return 0.0;

}
