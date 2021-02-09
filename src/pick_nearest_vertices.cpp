#include <pick_nearest_vertices.h>
#include <iostream>
#include <igl/Hit.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/unproject.h>

bool pick_nearest_vertices(std::vector<unsigned int> &verts, Eigen::Ref<const Eigen::Vector3d> win, 
                           Eigen::Ref<const Eigen::Matrix44f> view, Eigen::Ref<const Eigen::Matrix44f> proj, Eigen::Vector4f viewport,
                           Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double radius) {

    verts.clear();

    // Source, destination and direction in world
    Eigen::Vector3f start,dir;
    igl::Hit hit;

    //compute start and direction in the world to check for picked vertex
    //YOUR CODE HERE
    //
    // typedef Eigen::Matrix<typename Deriveds::Scalar,3,1> Vec3;
    // Vec3 win_s(win(0),win(1),0);
    // Vec3 win_d(win(0),win(1),1);
    // // Source, destination and direction in world
    // Vec3 d;

    Eigen::Vector3f win_0(win(0), win(1), win(2));
    Eigen::Vector3f win_1(win(0), win(1), 1.);


    igl::unproject(win_0,view,proj,viewport,start);
    igl::unproject(win_1,view,proj,viewport,dir);
    dir -= start;


    const auto & shoot_ray = [&V,&F](const Eigen::Vector3f& s, const Eigen::Vector3f& dir, igl::Hit & hit)->bool
    {
        std::vector<igl::Hit> hits;

        if(!igl::ray_mesh_intersect(s,dir,V,F,hits))
        {
            return false;
        }
        hit = hits[0];
        return true;
    };

    if(!shoot_ray(start,dir,hit))
    {
        return false;
    }

    //check if any of the hit vertices are within the selection radius
    //YOUR CODE HERE
    // Eigen::Vector3f hit_point = start + dir * hit.t;
    // for (int i=0; i<V.rows(); i++) {
    //     double sqr_dist = pow(V(i,0)-hit_point(0),2) + pow(V(i,1)-hit_point(1),2) + pow(V(i,2)-hit_point(2),2);
    //     if (sqr_dist <= pow(radius,2)) {
    //         verts.push_back(i);
    //     }
    // }

    Eigen::Vector3f bc;
    bc << 1.0-hit.u-hit.v, hit.u, hit.v;
    unsigned int fid = hit.id;

    long c;
    bc.maxCoeff(&c);
    unsigned int vid = F(fid,c);

    for(unsigned int qi = 0; qi < V.rows(); qi++) {
        if((V.row(qi) - V.row(vid)).norm() < radius) {

            verts.push_back(qi);
        }
    }

    return (verts.size() == 0 ? false : true);
}