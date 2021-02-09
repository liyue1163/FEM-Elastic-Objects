 
 #include <mass_matrix_linear_tetrahedron.h>

 void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
    Eigen::Matrix1212d M_partial;
    M_partial.setZero();
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            if (i == j) {
                M_partial.block<3, 3>(i*3, j*3) = 1.0/60 * Eigen::Matrix3d::Identity();
            }
            else {
                M_partial.block<3, 3>(i*3, j*3) = 1.0/120 * Eigen::Matrix3d::Identity();
            }
        }
    }

    M = 6.0 * density * volume * M_partial;
 }