#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T,
                                                   Eigen::Ref<const Eigen::MatrixXd> V_skin) {

    int n = V.rows();
    int M = T.rows();
    int L = V_skin.rows();

    N.resize(L, n);
    N.setZero();
    std::vector<Eigen::Triplet<double>> tripletList;

    for (int l = 0; l < L; l++) {
        const Eigen::Vector3d skin_x = V_skin.row(l);

        for (int m = 0; m < M; m++) {
            Eigen::Vector4d phi;
            const Eigen::RowVectorXi element = T.row(m);

            phi_linear_tetrahedron(phi, V, element, skin_x);
            // determine if point x on skin is inside current tetrahedron
            if (phi(0) >= 0 && phi(0) <= 1 &&
                phi(1) >= 0 && phi(1) <= 1 &&
                phi(2) >= 0 && phi(2) <= 1 &&
                phi(3) >= 0 && phi(3) <= 1 &&
                phi.sum() <= 1.05) {

                for (int i = 0; i < phi.size(); i++) {
                    tripletList.push_back({l, element(i), phi(i)});
                }

            }
        }
    }

    N.setFromTriplets(tripletList.begin(), tripletList.end());

}
