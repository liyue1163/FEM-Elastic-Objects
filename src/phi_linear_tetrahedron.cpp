#include <phi_linear_tetrahedron.h>

void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x) {
    Eigen::VectorXd X0 = V.row(element(0));
    Eigen::VectorXd X1 = V.row(element(1));
    Eigen::VectorXd X2 = V.row(element(2));
    Eigen::VectorXd X3 = V.row(element(3));

    Eigen::VectorXd dX = x - X0;

    Eigen::Matrix3d T;
    T.col(0) = (X1 - X0);
    T.col(1) = (X2 - X0);
    T.col(2) = (X3 - X0);

    Eigen::Vector3d phi_partial = T.inverse() * dX;
    phi(0) = 1.0 - phi_partial.sum();
    phi.tail<3>() = phi_partial;
}