#include <robot_init.h>
#include <cmath>

using namespace std;
using namespace ecn;

int main(int argc, char ** argv)
{
    // initialize robot class and get DOF's
    const auto robot{initRobot(argc, argv, 100)};
    const unsigned n{robot->getDofs()};

    // robot properties - max velocity and acceleration
    const auto vMax{robot->vMax()};
    const auto aMax{robot->aMax()};

    // main variables
    vpColVector q(n);               // joint position
    vpPoseVector p;                 // operational pose
    vpColVector qCommand(n);        // joint position setpoint
    vpColVector vCommand(n);        // joint velocity setpoint

    vpMatrix J;
    vpHomogeneousMatrix M;          // current pose
    vpHomogeneousMatrix M0, Md, Mi; // previous, final and current desired poses
    vpPoseVector pd;                // desired pose
    vpColVector v;                  // desired operational velocity
    vpColVector omega;              // desired operational angular velocity

    // TODO declare other variables if needed
    vpColVector q0(n), qf(n);       // joint position setpoint for initial and final poses
    double t0,tf,tfmax,alpha;       //misc doubles
    std::vector<double> tfv(n);     // vector to store tf max values for comparison
    vpMatrix R(6,6);                // declare empty R matrix for direct point-to-point control

    // main control loop
    while(robot->ok())
    {
        // current time
        const auto t{robot->time()};
        //std::cout<<(robot->mode()==ControlMode::VELOCITY_MANUAL);

        // update desired pose if has changed
        if(robot->newRef())
        {
            Md = robot->Md();
            //std::cout<<Md;
            M0 = robot->M0();
            pd.buildFrom(Md);
            t0 = t;
        }

        // get current joint positions
        q = robot->jointPosition();
        //cout << "Current joint position : " << q.t() << endl;

        // Direct Geometry for end-effector
        M = robot->fMe(q);  // matrix form
        p.buildFrom(M);     // translation + angle-axis form

        if(robot->mode() == ControlMode::POSITION_MANUAL)
        {
            // just check the Direct Geometric Model
            // TODO: fill the fMw function ** DONE
            robot->checkPose(M);
        }


        else if(robot->mode() == ControlMode::VELOCITY_MANUAL)
        {
            // follow a given operational velocity
            v = robot->guiVelocityScrew();

            // TODO: fill the fJw function **DONE
            // TODO: compute vCommand ------

            //using R, put rotation matrices from M on the diagonal of the 6x6 R matrix (declared at top)
            ecn::putAt(R,M.getRotationMatrix(),0,0);
            ecn::putAt(R,M.getRotationMatrix(),3,3);

            //qdot = J+*R*v
            vCommand = robot->fJe(q).pseudoInverse()*R*v;

            robot->setJointVelocity(vCommand);
        }


        else if(robot->mode() == ControlMode::DIRECT_P2P)
        {
            // find the Inverse Geometry to reach Md
            // TODO: fill the inverseGeometry function

            qf = robot->inverseGeometry(Md, q);
            robot->setJointPosition(qf);

        }

        else if(robot->mode() == ControlMode::POLYNOM_P2P)
        {
            // reach Md with polynomial joint trajectory
            // use q0 (initial position), qf (final), aMax and vMax

            // if reference has changed, compute new tf
            if(robot->newRef())
            {
                //set initial and final pose vectors
                q0 = robot->inverseGeometry(M0, q);
                qf = robot->inverseGeometry(Md, q);
                //set time at t0 to current robot time
                t0 = t;             

                //create vector of tfmax values for each joint
                for(unsigned int i =0; i<n; i++)
                    tfv[i] = max(abs(qf[i]-q0[i])*(3/(2*vMax[i])),sqrt((6*abs(qf[i]-q0[i]))/aMax[i]));

                //choose maximum tf from vector and apply to all joints
                tfmax = *max_element(tfv.begin(),tfv.end());
            }
            // TODO: compute qCommand from q0, qf, t, t0 and tf
            qCommand = q0 + (3*pow(min(((t-t0)/tfmax),1.0),2) - 2*pow(min(((t-t0)/tfmax),1.0),3))*(qf-q0);

            robot->setJointPosition(qCommand);
        }

        else if(robot->mode() == ControlMode::STRAIGHT_LINE_P2P)
        {
            // go from M0 to Md in 1 sec
            tf = 1.0;

            // TODO: compute qCommand from M0, Md, t, t0 and tf
            // use robot->intermediaryPose to build poses between M0 and Md

            alpha = fmin(t-t0,tf)/tf;

            qCommand = robot->inverseGeometry(robot->intermediaryPose(M0, Md, alpha),q);

            robot->setJointPosition(qCommand);
        }


        else if(robot->mode() == ControlMode::VELOCITY_P2P)
        {
            // go to Md using operational velocity

            v.resize(6,false);

            // TODO: compute joint velocity command

            pd.buildFrom(M.inverse()*Md);
            vpColVector theta_u = pd.getThetaUVector();
            auto t_trans = pd.getTranslationVector();

            ecn::putAt(v, M.getRotationMatrix()*t_trans,0);
            ecn::putAt(v, M.getRotationMatrix()*theta_u, 3);

            vCommand = robot->lambda()*robot->fJe(q).pseudoInverse()*v;

            robot->setJointVelocity(vCommand);
        }


    }

}
