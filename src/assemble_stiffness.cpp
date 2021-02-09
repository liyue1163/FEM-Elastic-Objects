#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0, 
                     double C, double D) {
    K.resize(q.rows(), q.rows());
    K.setZero();
    std::vector<Eigen::Triplet<double>> tripletList;

    // iterate through each tetrahedron
    for (int tet=0; tet<T.rows(); tet++) {

        Eigen::Matrix1212d K_tet;
        Eigen::RowVectorXi element = T.row(tet);
        double volume_i = v0(tet);
        d2V_linear_tetrahedron_dq2(K_tet, q, V, element, volume_i, C, D);

        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {

                Eigen::Matrix3d K_ij = -1.0 * K_tet.block<3,3>(i * 3,j * 3);
                int row_idx = element(i);
                int col_idx = element(j);

                for (int row_offset=0; row_offset<3; row_offset++) {
                    for (int col_offset=0; col_offset<3; col_offset++) {
                        tripletList.push_back({row_idx * 3 + row_offset, col_idx * 3 + col_offset, K_ij(row_offset, col_offset)});
                    }
                }
            }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());
        
    };
