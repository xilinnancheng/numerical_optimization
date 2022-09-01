#include "misc/visualizer.hpp"
#include "gcopter/cubic_curve.hpp"
#include "gcopter/cubic_spline.hpp"
#include "gcopter/path_smoother.hpp"

#include <ros/ros.h>
#include <ros/console.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include <random>

struct Config
{
    std::string targetTopic;
    double penaltyWeight;
    Eigen::Matrix3Xd circleObs;
    std::vector<Eigen::MatrixX4d> polygonObs;
    double pieceLength;
    double relCostTol;

    Config(const ros::NodeHandle &nh_priv)
    {
        std::vector<double> circleObsVec;
        std::vector<double> polygonObsVec;

        nh_priv.getParam("TargetTopic", targetTopic);
        nh_priv.getParam("PenaltyWeight", penaltyWeight);
        nh_priv.getParam("CircleObs", circleObsVec);
        nh_priv.getParam("PolygonObs", polygonObsVec);
        nh_priv.getParam("PieceLength", pieceLength);
        nh_priv.getParam("RelCostTol", relCostTol);

        circleObs = Eigen::Map<const Eigen::Matrix<double, 3, -1, Eigen::ColMajor>>(
            circleObsVec.data(), 3, circleObsVec.size() / 3);

        Eigen::Matrix<double, -1, 5, Eigen::RowMajor> polygonObs_tmp = Eigen::Map<const Eigen::Matrix<double, -1, 5, Eigen::RowMajor>>(
            polygonObsVec.data(), polygonObsVec.size() / 5, 5);

        if (polygonObs_tmp.rows() > 0)
        {
            double obs_id = polygonObs_tmp(0, 4);
            int start_row = 0;
            int end_row = 0;
            while (end_row < polygonObs_tmp.rows())
            {
                if (std::fabs(polygonObs_tmp(end_row, 4) - obs_id) > 1e-4)
                {
                    polygonObs.emplace_back(polygonObs_tmp.block(start_row, 0, end_row - start_row, 4));
                    start_row = end_row;
                    obs_id = polygonObs_tmp(end_row, 4);
                }
                end_row++;
            }

            if(end_row - start_row > 0){
                polygonObs.emplace_back(polygonObs_tmp.block(start_row, 0, end_row - start_row, 4));
            }
        }
    }
};

class CurveGen
{
private:
    Config config;
    ros::NodeHandle nh;
    ros::Subscriber targetSub;

    Visualizer visualizer;
    std::vector<Eigen::Vector2d> startGoal;

    CubicCurve curve;

public:
    CurveGen(ros::NodeHandle &nh_)
        : config(ros::NodeHandle("~")),
          nh(nh_),
          visualizer(nh)
    {
        targetSub = nh.subscribe(config.targetTopic, 1, &CurveGen::targetCallBack, this,
                                 ros::TransportHints().tcpNoDelay());
    }

    inline void vizObs()
    {
        // visualizer.visualizeDisks(config.circleObs);
        visualizer.visualizePolytope(config.polygonObs);
    }

    inline void plan()
    {
        if (startGoal.size() == 2)
        {
            const int N = (startGoal.back() - startGoal.front()).norm() / config.pieceLength;
            Eigen::Matrix2Xd innerPoints(2, N - 1);
            for (int i = 0; i < N - 1; ++i)
            {
                innerPoints.col(i) = (startGoal.back() - startGoal.front()) * (i + 1.0) / N + startGoal.front();
            }

            path_smoother::PathSmoother pathSmoother;
            pathSmoother.setup(startGoal.front(), startGoal.back(), N, config.circleObs, config.polygonObs, config.penaltyWeight);
            CubicCurve curve;

            if (std::isinf(pathSmoother.optimize(curve, innerPoints, config.relCostTol)))
            {
                return;
            }

            if (curve.getPieceNum() > 0)
            {
                visualizer.visualize(curve);
            }
        }
    }

    inline void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)
    {

        if (startGoal.size() >= 2)
        {
            startGoal.clear();
        }

        startGoal.emplace_back(msg->pose.position.x, msg->pose.position.y);

        plan();

        return;
    }
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "curve_gen_node");
    ros::NodeHandle nh_;

    CurveGen curveGen(nh_);

    ros::Duration(2.0).sleep();

    curveGen.vizObs();

    ros::spin();

    return 0;
}
