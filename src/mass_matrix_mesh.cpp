#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {

    M.resize(qdot.rows(), qdot.rows());
    std::vector<Eigen::Triplet<double>> tripletList;

    for (int tet=0; tet<T.rows(); tet++) {
        Eigen::Matrix1212d M_j;
        Eigen::RowVectorXi element_j = T.row(tet);
        double volume_j = v0(tet);
        mass_matrix_linear_tetrahedron(M_j, qdot, element_j, density, volume_j);

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int row_idx = element_j(i);
                int col_idx = element_j(j);

                for (int k=0; k<3; k++) {
                    tripletList.push_back({row_idx * 3 + k, col_idx * 3 + k, M_j(i * 3 + k, j * 3 + k)});
                }
            }
        }
    }
    M.setFromTriplets(tripletList.begin(), tripletList.end());

}