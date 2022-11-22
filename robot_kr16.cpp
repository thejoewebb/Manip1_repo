#include <robot_kr16.h>
#include <trig_solvers.h>

// Model of Kuka KR16 robot

auto Re = 0.158;
const auto r1{0.675};
const auto a2{0.26};
const auto a3{0.68};
const auto a4{0.035};
const auto r4{0.67};

// Any end-effector to wrist constant transform
void ecn::RobotKr16::init_wMe()
{

    // Generated end-effector code
        wMe[0][0] = 1.;
        wMe[0][1] = 0;
        wMe[0][2] = 0;
        wMe[0][3] = 0;
        wMe[1][0] = 0;
        wMe[1][1] = -1.;
        wMe[1][2] = 0;
        wMe[1][3] = 0;
        wMe[2][0] = 0;
        wMe[2][1] = 0;
        wMe[2][2] = -1.;
        wMe[2][3] = -Re;
        wMe[3][0] = 0;
        wMe[3][1] = 0;
        wMe[3][2] = 0;
        wMe[3][3] = 1.;
        // End of end-effector code
}

// Direct Geometry fixed to wrist frame
vpHomogeneousMatrix ecn::RobotKr16::fMw(const vpColVector &q) const
{
    vpHomogeneousMatrix M;

    // Generated pose code
        const double c1 = cos(q[0]);
        const double c2 = cos(q[1]);
        const double c4 = cos(q[3]);
        const double c5 = cos(q[4]);
        const double c6 = cos(q[5]);
        const double c23 = cos(q[1]+q[2]);
        const double s1 = sin(q[0]);
        const double s2 = sin(q[1]);
        const double s4 = sin(q[3]);
        const double s5 = sin(q[4]);
        const double s6 = sin(q[5]);
        const double s23 = sin(q[1]+q[2]);
        M[0][0] = (-(s1*s4 + s23*c1*c4)*c5 - s5*c1*c23)*c6 - (s1*c4 - s4*s23*c1)*s6;
        M[0][1] = -(-(s1*s4 + s23*c1*c4)*c5 - s5*c1*c23)*s6 - (s1*c4 - s4*s23*c1)*c6;
        M[0][2] = (s1*s4 + s23*c1*c4)*s5 - c1*c5*c23;
        M[0][3] = (a2 + a3*c2 - a4*s23 + r4*c23)*c1;
        M[1][0] = ((s1*s23*c4 - s4*c1)*c5 + s1*s5*c23)*c6 - (s1*s4*s23 + c1*c4)*s6;
        M[1][1] = -((s1*s23*c4 - s4*c1)*c5 + s1*s5*c23)*s6 - (s1*s4*s23 + c1*c4)*c6;
        M[1][2] = -(s1*s23*c4 - s4*c1)*s5 + s1*c5*c23;
        M[1][3] = (-a2 - a3*c2 + a4*s23 - r4*c23)*s1;
        M[2][0] = (s5*s23 - c4*c5*c23)*c6 + s4*s6*c23;
        M[2][1] = -(s5*s23 - c4*c5*c23)*s6 + s4*c6*c23;
        M[2][2] = s5*c4*c23 + s23*c5;
        M[2][3] = -a3*s2 - a4*c23 + r1 - r4*s23;
        M[3][0] = 0;
        M[3][1] = 0;
        M[3][2] = 0;
        M[3][3] = 1.;
    // End of pose code

    return M;
}


// Inverse Geometry
vpColVector ecn::RobotKr16::inverseGeometry(const vpHomogeneousMatrix &fMe_des, const vpColVector &q0) const
{
    // build corresponding oMw and explode into 12 elements
    const auto [xx,xy,xz,yx,yy,yz,zx,zy,zz,tx,ty,tz] = explodeMatrix(fMe_des);

    double alpha;
    vpRotationMatrix R_03;
    auto M_06 = fMe_des*wMe.inverse();

    // TODO add candidates

    for(const double q1: solveType2(tx,ty,0)){

        const auto c1{cos(q1)};
        const auto s1{sin(q1)};

        if(c1>s1)
            alpha = tx/c1;
        else
            alpha = -ty/s1;

        for(const auto [q2,q23]: solveType7(0, -a3, (r1-tz), (a2-alpha), a4, r4)){

            const auto q3{q23-q2};
            const auto c23{cos(q23)};
            const auto s23{sin(q23)};

            //build rotation matrix from base to joint 3
            R_03[0][0] = -s23*c1;
            R_03[0][1] = -c1*c23;
            R_03[0][2] = s1;
            R_03[1][0] = s1*s23;
            R_03[1][1] = s1*c23;
            R_03[1][2] = c1;
            R_03[2][0] = -c23;
            R_03[2][1] = s23;
            R_03[2][2] = 0;

            //compute rotation matrix from wrist to end effector frame
            auto R_36 = R_03.t()*M_06.getRotationMatrix();

            //explodematrix again
            auto zx2 = R_36[0][2];
            auto zy2 = R_36[1][2];
            auto zz2 = R_36[2][2];
            auto xx2 = R_36[0][0];
            auto xy2 = R_36[1][0];
            auto xz2 = R_36[2][0];
            auto yy2 = R_36[1][1];

            if(isNull(zx2||zz2)){
                const auto q5{0};
                const auto q4{q0[3]};
                const auto q46 = atan2(-xz2,-xx2);
                const auto q6 = q46 - q4;

                addCandidate({q1,q2,q3,q4,q5,q6});

            }else{
                for(const double q5: solveType2(0,1,zy2)){

                    auto s5 = sin(q5);

                    for(const double q4: solveType3(0,-s5,zx2,s5,0,zz2)){

                        for(const double q6: solveType2(0,s5,xy2)){

                            addCandidate({q1,q2,q3,q4,q5,q6});
                        }
                    }
                }
            }
        }
    }

    return bestCandidate(q0);
}


vpMatrix ecn::RobotKr16::fJw(const vpColVector &q) const
{
    vpMatrix J(6, dofs);

    // Generated Jacobian code
        const double c1 = cos(q[0]);
        const double c2 = cos(q[1]);
        const double c4 = cos(q[3]);
        const double c5 = cos(q[4]);
        const double c23 = cos(q[1]+q[2]);
        const double s1 = sin(q[0]);
        const double s2 = sin(q[1]);
        const double s4 = sin(q[3]);
        const double s5 = sin(q[4]);
        const double s23 = sin(q[1]+q[2]);
        J[0][0] = (0.035*s23 - 0.68*c2 - 0.67*c23 - 0.26)*s1;
        J[0][1] = -(0.68*s2 + 0.67*s23 + 0.035*c23)*c1;
        J[0][2] = -(0.67*s23 + 0.035*c23)*c1;
        //J[0][3] = 0;
        //J[0][4] = 0;
        //J[0][5] = 0;
        J[1][0] = -(-0.035*s23 + 0.68*c2 + 0.67*c23 + 0.26)*c1;
        J[1][1] = (0.68*s2 + 0.67*s23 + 0.035*c23)*s1;
        J[1][2] = (0.67*s23 + 0.035*c23)*s1;
        //J[1][3] = 0;
        //J[1][4] = 0;
        //J[1][5] = 0;
        //J[2][0] = 0;
        J[2][1] = 0.035*s23 - 0.68*c2 - 0.67*c23;
        J[2][2] = 0.035*s23 - 0.67*c23;
        //J[2][3] = 0;
        //J[2][4] = 0;
        //J[2][5] = 0;
        //J[3][0] = 0;
        J[3][1] = s1;
        J[3][2] = s1;
        J[3][3] = -c1*c23;
        J[3][4] = s1*c4 - s4*s23*c1;
        J[3][5] = (s1*s4 + s23*c1*c4)*s5 - c1*c5*c23;
        //J[4][0] = 0;
        J[4][1] = c1;
        J[4][2] = c1;
        J[4][3] = s1*c23;
        J[4][4] = s1*s4*s23 + c1*c4;
        J[4][5] = -(s1*s23*c4 - s4*c1)*s5 + s1*c5*c23;
        J[5][0] = -1.;
        //J[5][1] = 0;
        //J[5][2] = 0;
        J[5][3] = s23;
        J[5][4] = -s4*c23;
        J[5][5] = s5*c4*c23 + s23*c5;
        // End of Jacobian code

    return J;
}
