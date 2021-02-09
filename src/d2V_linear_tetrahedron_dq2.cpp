#include <d2V_linear_tetrahedron_dq2.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <d2psi_neo_hookean_dq2.h>
#include <quadrature_single_point.h>

void d2V_linear_tetrahedron_dq2(Eigen::Matrix1212d &H, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

   auto neohookean_linear_tet = [&](Eigen::Matrix1212d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        
      //Code to compute non-integrated hessian matrix goes here
       Eigen::MatrixXd x;
       x.resize(3, 4);
       for (int i = 0; i < 4; i++) {
           x.col(i) = q.segment(element(i) * 3, 3);
       }

       Eigen::Matrix43d dphi;
       dphi_linear_tetrahedron_dX(dphi, V, element, X);

       Eigen::Matrix3d F;
       F = x * dphi;

       Eigen::Matrix99d hessian_F;
       d2psi_neo_hookean_dF2(hessian_F, F, C, D);

       Eigen::MatrixXd B_j;
       B_j.resize(9, 12);
       B_j.setZero();

       // 1st 9x12 block
       B_j(0, 0) = dphi(0, 0);
       B_j(1, 0) = dphi(0, 1);
       B_j(2, 0) = dphi(0, 2);

       B_j(3, 1) = dphi(0, 0);
       B_j(4, 1) = dphi(0, 1);
       B_j(5, 1) = dphi(0, 2);

       B_j(6, 2) = dphi(0, 0);
       B_j(7, 2) = dphi(0, 1);
       B_j(8, 2) = dphi(0, 2);

       // 2nd 9x12 block
       B_j(0, 3) = dphi(1, 0);
       B_j(1, 3) = dphi(1, 1);
       B_j(2, 3) = dphi(1, 2);

       B_j(3, 4) = dphi(1, 0);
       B_j(4, 4) = dphi(1, 1);
       B_j(5, 4) = dphi(1, 2);

       B_j(6, 5) = dphi(1, 0);
       B_j(7, 5) = dphi(1, 1);
       B_j(8, 5) = dphi(1, 2);

       // 3rd 9x12 block
       B_j(0, 6) = dphi(2, 0);
       B_j(1, 6) = dphi(2, 1);
       B_j(2, 6) = dphi(2, 2);

       B_j(3, 7) = dphi(2, 0);
       B_j(4, 7) = dphi(2, 1);
       B_j(5, 7) = dphi(2, 2);

       B_j(6, 8) = dphi(2, 0);
       B_j(7, 8) = dphi(2, 1);
       B_j(8, 8) = dphi(2, 2);

       // 4th 9x12 block
       B_j(0, 9) = dphi(3, 0);
       B_j(1, 9) = dphi(3, 1);
       B_j(2, 9) = dphi(3, 2);

       B_j(3, 10) = dphi(3, 0);
       B_j(4, 10) = dphi(3, 1);
       B_j(5, 10) = dphi(3, 2);

       B_j(6, 11) = dphi(3, 0);
       B_j(7, 11) = dphi(3, 1);
       B_j(8, 11) = dphi(3, 2);

       dV = B_j.transpose() * hessian_F * B_j;

    };

    //integrate the non-integrated hessian across the tetrahedral element
    quadrature_single_point(H, q, element, volume, neohookean_linear_tet);  
    

    //DO NOT REMOVE THIS CODE This code ensures that the hessian matrix is symmetric postive definite by projecting all
    //negative eigenvalues to small, postive values.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix1212d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 12; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }

    H = Evec * DiagEval * Evec.transpose();

}
