#include <robot_turret.h>
#include <trig_solvers.h>

// Model of Turret robot

auto r1{0.5};
auto r2{0.1};
auto re{0};

// Any constant transform at base or end-effector
void ecn::RobotTurret::init_wMe()
{

    // Generated end-effector code
        wMe[0][0] = 1.;
        wMe[0][1] = 0;
        wMe[0][2] = 0;
        wMe[0][3] = 0;
        wMe[1][0] = 0;
        wMe[1][1] = 1.;
        wMe[1][2] = 0;
        wMe[1][3] = 0;
        wMe[2][0] = 0;
        wMe[2][1] = 0;
        wMe[2][2] = 1.;
        wMe[2][3] = re;
        wMe[3][0] = 0;
        wMe[3][1] = 0;
        wMe[3][2] = 0;
        wMe[3][3] = 1.;
        // End of end-effector code

}


// Direct Geometry
vpHomogeneousMatrix ecn::RobotTurret::fMw(const vpColVector &q) const
{

    vpHomogeneousMatrix M;

    // Generated pose code
        const double c1 = cos(q[0]);
        const double c2 = cos(q[1]);
        const double s1 = sin(q[0]);
        const double s2 = sin(q[1]);
       M[0][0] = c1*c2;
       M[0][1] = -s1;
       M[0][2] = -s2*c1;
       M[0][3] = (-q[2] - r2)*s2*c1;
       M[1][0] = s1*c2;
       M[1][1] = c1;
       M[1][2] = -s1*s2;
       M[1][3] = (-q[2] - r2)*s1*s2;
       M[2][0] = s2;
       M[2][2] = c2;
       M[2][3] = r1 + (q[2] + r2)*c2;
       M[3][3] = 1.;
       // End of pose code

    return M;
}


// Inverse Geometry
vpColVector ecn::RobotTurret::inverseGeometry(const vpHomogeneousMatrix &fMe_des, const vpColVector &q0) const
{
    // build corresponding oMw and explode into 12 elements
    const auto [xx,xy,xz,yx,yy,yz,zx,zy,zz,tx,ty,tz] = explodeMatrix(fMe_des);
    //double q1, q2, q3;

    // TODO add candidates

    for(const auto q1: solveType3(-1,0,yx,0,1,yy)){

        const auto c1{cos(q1)};
        const auto s1{sin(q1)};

        for(const auto q2: solveType3(1,0,xz,0,1,zz)){

            const auto c2{cos(q2)};
            const auto s2{sin(q2)};

            if(!isNull(c2)){

                const auto q3 = ((tz-r1)/c2)-r2;

                addCandidate({q1,q2,q3});

            }else if(!isNull(s1)){

                const auto q3 = -(ty/(s1*s2)+r2);

                addCandidate({q1,q2,q3});

            }else{
                const auto q3 = (-tx/s2*c1)-r2;

                addCandidate({q1,q2,q3});
            }
        }
    }

return bestCandidate(q0);
}

// Wrist Jacobian
vpMatrix ecn::RobotTurret::fJw(const vpColVector &q) const
{
    vpMatrix J(6, dofs);

    // Generated Jacobian code
        const double c1 = cos(q[0]);
        const double c2 = cos(q[1]);
        const double s1 = sin(q[0]);
        const double s2 = sin(q[1]);
        J[0][0] = (q[2] + r2)*s1*s2;
        J[0][1] = (-q[2] - r2)*c1*c2;
        J[0][2] = -s2*c1;
        J[1][0] = (-q[2] - r2)*s2*c1;
        J[1][1] = (-q[2] - r2)*s1*c2;
        J[1][2] = -s1*s2;
        J[2][1] = -(q[2] + r2)*s2;
        J[2][2] = c2;
        J[3][1] = s1;
        J[4][1] = -c1;
        J[5][0] = 1.;

        // End of Jacobian code

    return J;
}
