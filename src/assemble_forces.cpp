#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D) {
    f = Eigen::VectorXd::Zero(q.rows());

    // iterate through each tetrahedron
    for (int i=0; i<T.rows(); i++) {

        Eigen::Vector12d dV;
        Eigen::RowVectorXi element = T.row(i);
        double volume_i = v0(i);
        dV_linear_tetrahedron_dq(dV, q, V, element, volume_i, C, D);

        int tet_idx1 = element(0);
        int tet_idx2 = element(1);
        int tet_idx3 = element(2);
        int tet_idx4 = element(3);

        // directly update relevant elements from f,
        // bunny moves faster and more smoothly
        f[3*tet_idx1] -= dV[0];
        f[3*tet_idx1+1] -= dV[1];
        f[3*tet_idx1+2] -= dV[2];

        f[3*tet_idx2] -= dV[3];
        f[3*tet_idx2+1] -= dV[4];
        f[3*tet_idx2+2] -= dV[5];

        f[3*tet_idx3] -= dV[6];
        f[3*tet_idx3+1] -= dV[7];
        f[3*tet_idx3+2] -= dV[8];

        f[3*tet_idx4] -= dV[9];
        f[3*tet_idx4+1] -= dV[10];
        f[3*tet_idx4+2] -= dV[11];

    }
};