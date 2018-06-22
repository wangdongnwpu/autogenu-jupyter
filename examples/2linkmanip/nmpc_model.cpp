#include "nmpc_model.hpp"



void nmpc_model::phix(const double t, const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> lmd)
{
	lmd[0] = (x[0] - xf[0]) * sf[0];
	lmd[1] = (x[1] - xf[1]) * sf[1];
	lmd[2] = (x[2] - xf[2]) * sf[2];
	lmd[3] = (x[3] - xf[3]) * sf[3];
}

void nmpc_model::statefunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> f)
{
	f[0] = x[2];
	f[1] = x[3];
	f[2] = (l1 * m2 * m2 * (0.2e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) + g) * sin(x[0]) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - l1 * m2 * (m2 * (0.2e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - x[2] * l1 * (x[2] + x[3])) * d2 * sin(x[0] + x[1]) + x[2] * cos(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 + x[2] * cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + (d2 * d2 * (x[2] + x[3]) * m2 + J2 * x[2] + J2 * x[3] - u[1]) * cos(x[0]) - x[2] * sin(x[0]) * (d2 * d2 * m2 + J2) * (x[2] + x[3])) * d2 * cos(x[0] + x[1]) + l1 * m2 * (-x[2] * cos(x[2]) * sin(x[0]) * d2 * l1 * m2 * cos(x[2] + x[3]) - x[2] * sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 - x[2] * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * cos(x[0]) + sin(x[0]) * (-d2 * d2 * (x[2] + x[3]) * m2 - J2 * x[2] - J2 * x[3] + u[1])) * d2 * sin(x[0] + x[1]) + cos(x[2]) * d2 * l1 * m2 * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * cos(x[2] + x[3]) + sin(x[2]) * d2 * l1 * m2 * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * sin(x[2] + x[3]) - x[2] * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * (x[2] + x[3]) * cos(x[0]) + ((d1 * d2 * d2 * m1 + J2 * l1) * m2 + J2 * d1 * m1) * g * sin(x[0]) + (d2 * d2 * m2 + J2) * (l1 * l1 * m2 * x[2] + (d1 * d1 * m1 + J1) * x[2] + u[1] - u[0])) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1));
	f[3] = (-(0.4e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) + g) * l1 * m2 * m2 * sin(x[0]) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - l1 * m2 * (-m2 * (0.4e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - 0.2e1 * x[2] * l1 * (x[2] + x[3])) * d2 * sin(x[0] + x[1]) + x[3] * cos(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 + x[3] * cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + (g * (d1 * m1 + l1 * m2) * sin(x[0]) + ((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3] + 0.2e1 * u[1] - u[0]) * cos(x[0]) + ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * x[2] * sin(x[0])) * d2 * cos(x[0] + x[1]) + m2 * (-x[3] * cos(x[2]) * sin(x[0]) * d2 * l1 * l1 * m2 * cos(x[2] + x[3]) - x[3] * sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * l1 * m2 + g * l1 * (d1 * m1 + l1 * m2) * pow(cos(x[0]), 0.2e1) + l1 * ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * x[2] * cos(x[0]) - 0.2e1 * l1 * (((-x[2] / 0.2e1 - x[3] / 0.2e1) * d2 * d2 + l1 * l1 * x[2] / 0.2e1) * m2 + (d1 * d1 * m1 / 0.2e1 + J1 / 0.2e1 - J2 / 0.2e1) * x[2] - J2 * x[3] / 0.2e1 + u[1] - u[0] / 0.2e1) * sin(x[0]) + g * (d1 * d1 * m1 - d1 * l1 * m1 + J1)) * d2 * sin(x[0] + x[1]) + l1 * m2 * (((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3]) * cos(x[2]) * d2 * cos(x[2] + x[3]) + l1 * sin(x[2]) * m2 * (((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3]) * d2 * sin(x[2] + x[3]) + 0.2e1 * x[2] * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * (x[2] + x[3]) * cos(x[0]) - ((d1 * d2 * d2 * m1 + J2 * l1) * m2 + J2 * d1 * m1) * g * sin(x[0]) + d2 * d2 * l1 * l1 * m2 * m2 * x[3] + (((d1 * d1 * m1 + J1) * x[3] - u[1] + u[0]) * d2 * d2 - l1 * l1 * (-J2 * x[3] + u[1])) * m2 + J2 * (d1 * d1 * m1 + J1) * x[3] + (-u[1] + u[0]) * J2 - u[1] * (d1 * d1 * m1 + J1)) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1));

}


void nmpc_model::hxfunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hx)
{
	hx[0] = -( -(x[0] - xf[0]) * q[0] - lmd[2] * (-0.2e1 * l1 * l1 * m2 * m2 * x[2] * (x[2] + x[3]) * pow(sin(x[0]), 0.2e1) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + l1 * m2 * m2 * (0.2e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) + g) * cos(x[0]) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - 0.2e1 * l1 * m2 * m2 * (0.2e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) + g) * sin(x[0]) * d2 * d2 * cos(x[0] + x[1]) * sin(x[0] + x[1]) - l1 * m2 * (m2 * (-0.4e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) * sin(x[0]) - g * sin(x[0])) * d2 * sin(x[0] + x[1]) + m2 * (0.2e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - x[2] * l1 * (x[2] + x[3])) * d2 * cos(x[0] + x[1]) - x[2] * cos(x[2]) * sin(x[0]) * d2 * l1 * m2 * cos(x[2] + x[3]) - x[2] * sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 - (d2 * d2 * (x[2] + x[3]) * m2 + J2 * x[2] + J2 * x[3] - u[1]) * sin(x[0]) - x[2] * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * cos(x[0])) * d2 * cos(x[0] + x[1]) + l1 * m2 * (m2 * (0.2e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - x[2] * l1 * (x[2] + x[3])) * d2 * sin(x[0] + x[1]) + x[2] * cos(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 + x[2] * cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + (d2 * d2 * (x[2] + x[3]) * m2 + J2 * x[2] + J2 * x[3] - u[1]) * cos(x[0]) - x[2] * sin(x[0]) * (d2 * d2 * m2 + J2) * (x[2] + x[3])) * d2 * sin(x[0] + x[1]) + l1 * m2 * (-x[2] * cos(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 - x[2] * cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + x[2] * sin(x[0]) * (d2 * d2 * m2 + J2) * (x[2] + x[3]) + cos(x[0]) * (-d2 * d2 * (x[2] + x[3]) * m2 - J2 * x[2] - J2 * x[3] + u[1])) * d2 * sin(x[0] + x[1]) + l1 * m2 * (-x[2] * cos(x[2]) * sin(x[0]) * d2 * l1 * m2 * cos(x[2] + x[3]) - x[2] * sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 - x[2] * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * cos(x[0]) + sin(x[0]) * (-d2 * d2 * (x[2] + x[3]) * m2 - J2 * x[2] - J2 * x[3] + u[1])) * d2 * cos(x[0] + x[1]) - x[2] * pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 * (x[2] + x[3]) + x[2] * pow(sin(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 * (x[2] + x[3]) + ((d1 * d2 * d2 * m1 + J2 * l1) * m2 + J2 * d1 * m1) * g * cos(x[0])) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1)) + lmd[2] * (l1 * m2 * m2 * (0.2e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) + g) * sin(x[0]) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - l1 * m2 * (m2 * (0.2e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - x[2] * l1 * (x[2] + x[3])) * d2 * sin(x[0] + x[1]) + x[2] * cos(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 + x[2] * cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + (d2 * d2 * (x[2] + x[3]) * m2 + J2 * x[2] + J2 * x[3] - u[1]) * cos(x[0]) - x[2] * sin(x[0]) * (d2 * d2 * m2 + J2) * (x[2] + x[3])) * d2 * cos(x[0] + x[1]) + l1 * m2 * (-x[2] * cos(x[2]) * sin(x[0]) * d2 * l1 * m2 * cos(x[2] + x[3]) - x[2] * sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 - x[2] * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * cos(x[0]) + sin(x[0]) * (-d2 * d2 * (x[2] + x[3]) * m2 - J2 * x[2] - J2 * x[3] + u[1])) * d2 * sin(x[0] + x[1]) + cos(x[2]) * d2 * l1 * m2 * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * cos(x[2] + x[3]) + sin(x[2]) * d2 * l1 * m2 * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * sin(x[2] + x[3]) - x[2] * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * (x[2] + x[3]) * cos(x[0]) + ((d1 * d2 * d2 * m1 + J2 * l1) * m2 + J2 * d1 * m1) * g * sin(x[0]) + (d2 * d2 * m2 + J2) * (l1 * l1 * m2 * x[2] + (d1 * d1 * m1 + J1) * x[2] + u[1] - u[0])) * pow((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1), -0.2e1) * (-0.2e1 * cos(x[0]) * sin(x[0]) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - 0.2e1 * (0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * cos(x[0] + x[1]) * sin(x[0] + x[1]) - 0.2e1 * pow(sin(x[0]), 0.2e1) * cos(x[0] + x[1]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - 0.2e1 * cos(x[0]) * pow(sin(x[0] + x[1]), 0.2e1) * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 + 0.2e1 * pow(cos(x[0]), 0.2e1) * cos(x[0] + x[1]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 + 0.2e1 * cos(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * sin(x[0])) - lmd[3] * (0.4e1 * l1 * l1 * m2 * m2 * x[2] * (x[2] + x[3]) * pow(sin(x[0]), 0.2e1) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - (0.4e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) + g) * l1 * m2 * m2 * cos(x[0]) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * (0.4e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) + g) * l1 * m2 * m2 * sin(x[0]) * d2 * d2 * cos(x[0] + x[1]) * sin(x[0] + x[1]) - l1 * m2 * (-m2 * (-0.8e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) * sin(x[0]) - g * sin(x[0])) * d2 * sin(x[0] + x[1]) - m2 * (0.4e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - 0.2e1 * x[2] * l1 * (x[2] + x[3])) * d2 * cos(x[0] + x[1]) - x[3] * sin(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 - x[3] * sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + g * (d1 * m1 + l1 * m2) * pow(cos(x[0]), 0.2e1) - (g * (d1 * m1 + l1 * m2) * sin(x[0]) + ((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3] + 0.2e1 * u[1] - u[0]) * sin(x[0]) + ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * x[2] * cos(x[0])) * d2 * cos(x[0] + x[1]) + l1 * m2 * (-m2 * (0.4e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - 0.2e1 * x[2] * l1 * (x[2] + x[3])) * d2 * sin(x[0] + x[1]) + x[3] * cos(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 + x[3] * cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + (g * (d1 * m1 + l1 * m2) * sin(x[0]) + ((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3] + 0.2e1 * u[1] - u[0]) * cos(x[0]) + ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * x[2] * sin(x[0])) * d2 * sin(x[0] + x[1]) + m2 * (-x[3] * cos(x[2]) * cos(x[0]) * d2 * l1 * l1 * m2 * cos(x[2] + x[3]) - x[3] * cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * l1 * m2 - 0.2e1 * g * l1 * (d1 * m1 + l1 * m2) * cos(x[0]) * sin(x[0]) - l1 * ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * x[2] * sin(x[0]) - 0.2e1 * l1 * (((-x[2] / 0.2e1 - x[3] / 0.2e1) * d2 * d2 + l1 * l1 * x[2] / 0.2e1) * m2 + (d1 * d1 * m1 / 0.2e1 + J1 / 0.2e1 - J2 / 0.2e1) * x[2] - J2 * x[3] / 0.2e1 + u[1] - u[0] / 0.2e1) * cos(x[0])) * d2 * sin(x[0] + x[1]) + m2 * (-x[3] * cos(x[2]) * sin(x[0]) * d2 * l1 * l1 * m2 * cos(x[2] + x[3]) - x[3] * sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * l1 * m2 + g * l1 * (d1 * m1 + l1 * m2) * pow(cos(x[0]), 0.2e1) + l1 * ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * x[2] * cos(x[0]) - 0.2e1 * l1 * (((-x[2] / 0.2e1 - x[3] / 0.2e1) * d2 * d2 + l1 * l1 * x[2] / 0.2e1) * m2 + (d1 * d1 * m1 / 0.2e1 + J1 / 0.2e1 - J2 / 0.2e1) * x[2] - J2 * x[3] / 0.2e1 + u[1] - u[0] / 0.2e1) * sin(x[0]) + g * (d1 * d1 * m1 - d1 * l1 * m1 + J1)) * d2 * cos(x[0] + x[1]) + 0.2e1 * x[2] * pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 * (x[2] + x[3]) - 0.2e1 * x[2] * pow(sin(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 * (x[2] + x[3]) - ((d1 * d2 * d2 * m1 + J2 * l1) * m2 + J2 * d1 * m1) * g * cos(x[0])) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1)) + lmd[3] * (-(0.4e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) + g) * l1 * m2 * m2 * sin(x[0]) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - l1 * m2 * (-m2 * (0.4e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - 0.2e1 * x[2] * l1 * (x[2] + x[3])) * d2 * sin(x[0] + x[1]) + x[3] * cos(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 + x[3] * cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + (g * (d1 * m1 + l1 * m2) * sin(x[0]) + ((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3] + 0.2e1 * u[1] - u[0]) * cos(x[0]) + ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * x[2] * sin(x[0])) * d2 * cos(x[0] + x[1]) + m2 * (-x[3] * cos(x[2]) * sin(x[0]) * d2 * l1 * l1 * m2 * cos(x[2] + x[3]) - x[3] * sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * l1 * m2 + g * l1 * (d1 * m1 + l1 * m2) * pow(cos(x[0]), 0.2e1) + l1 * ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * x[2] * cos(x[0]) - 0.2e1 * l1 * (((-x[2] / 0.2e1 - x[3] / 0.2e1) * d2 * d2 + l1 * l1 * x[2] / 0.2e1) * m2 + (d1 * d1 * m1 / 0.2e1 + J1 / 0.2e1 - J2 / 0.2e1) * x[2] - J2 * x[3] / 0.2e1 + u[1] - u[0] / 0.2e1) * sin(x[0]) + g * (d1 * d1 * m1 - d1 * l1 * m1 + J1)) * d2 * sin(x[0] + x[1]) + l1 * m2 * (((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3]) * cos(x[2]) * d2 * cos(x[2] + x[3]) + l1 * sin(x[2]) * m2 * (((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3]) * d2 * sin(x[2] + x[3]) + 0.2e1 * x[2] * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * (x[2] + x[3]) * cos(x[0]) - ((d1 * d2 * d2 * m1 + J2 * l1) * m2 + J2 * d1 * m1) * g * sin(x[0]) + d2 * d2 * l1 * l1 * m2 * m2 * x[3] + (((d1 * d1 * m1 + J1) * x[3] - u[1] + u[0]) * d2 * d2 - l1 * l1 * (-J2 * x[3] + u[1])) * m2 + J2 * (d1 * d1 * m1 + J1) * x[3] + (-u[1] + u[0]) * J2 - u[1] * (d1 * d1 * m1 + J1)) * pow((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1), -0.2e1) * (-0.2e1 * cos(x[0]) * sin(x[0]) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - 0.2e1 * (0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * cos(x[0] + x[1]) * sin(x[0] + x[1]) - 0.2e1 * pow(sin(x[0]), 0.2e1) * cos(x[0] + x[1]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - 0.2e1 * cos(x[0]) * pow(sin(x[0] + x[1]), 0.2e1) * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 + 0.2e1 * pow(cos(x[0]), 0.2e1) * cos(x[0] + x[1]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 + 0.2e1 * cos(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * sin(x[0])) );
	hx[1] = -( -(x[1] - xf[1]) * q[1] - lmd[2] * (-0.2e1 * l1 * m2 * m2 * (0.2e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) + g) * sin(x[0]) * d2 * d2 * cos(x[0] + x[1]) * sin(x[0] + x[1]) - l1 * m2 * m2 * (0.2e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - x[2] * l1 * (x[2] + x[3])) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + l1 * m2 * (m2 * (0.2e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - x[2] * l1 * (x[2] + x[3])) * d2 * sin(x[0] + x[1]) + x[2] * cos(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 + x[2] * cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + (d2 * d2 * (x[2] + x[3]) * m2 + J2 * x[2] + J2 * x[3] - u[1]) * cos(x[0]) - x[2] * sin(x[0]) * (d2 * d2 * m2 + J2) * (x[2] + x[3])) * d2 * sin(x[0] + x[1]) + l1 * m2 * (-x[2] * cos(x[2]) * sin(x[0]) * d2 * l1 * m2 * cos(x[2] + x[3]) - x[2] * sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 - x[2] * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * cos(x[0]) + sin(x[0]) * (-d2 * d2 * (x[2] + x[3]) * m2 - J2 * x[2] - J2 * x[3] + u[1])) * d2 * cos(x[0] + x[1])) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1)) + lmd[2] * (l1 * m2 * m2 * (0.2e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) + g) * sin(x[0]) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - l1 * m2 * (m2 * (0.2e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - x[2] * l1 * (x[2] + x[3])) * d2 * sin(x[0] + x[1]) + x[2] * cos(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 + x[2] * cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + (d2 * d2 * (x[2] + x[3]) * m2 + J2 * x[2] + J2 * x[3] - u[1]) * cos(x[0]) - x[2] * sin(x[0]) * (d2 * d2 * m2 + J2) * (x[2] + x[3])) * d2 * cos(x[0] + x[1]) + l1 * m2 * (-x[2] * cos(x[2]) * sin(x[0]) * d2 * l1 * m2 * cos(x[2] + x[3]) - x[2] * sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 - x[2] * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * cos(x[0]) + sin(x[0]) * (-d2 * d2 * (x[2] + x[3]) * m2 - J2 * x[2] - J2 * x[3] + u[1])) * d2 * sin(x[0] + x[1]) + cos(x[2]) * d2 * l1 * m2 * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * cos(x[2] + x[3]) + sin(x[2]) * d2 * l1 * m2 * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * sin(x[2] + x[3]) - x[2] * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * (x[2] + x[3]) * cos(x[0]) + ((d1 * d2 * d2 * m1 + J2 * l1) * m2 + J2 * d1 * m1) * g * sin(x[0]) + (d2 * d2 * m2 + J2) * (l1 * l1 * m2 * x[2] + (d1 * d1 * m1 + J1) * x[2] + u[1] - u[0])) * pow((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1), -0.2e1) * (-0.2e1 * (0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * cos(x[0] + x[1]) * sin(x[0] + x[1]) - 0.2e1 * cos(x[0]) * pow(sin(x[0] + x[1]), 0.2e1) * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 + 0.2e1 * cos(x[0]) * sin(x[0]) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1)) - lmd[3] * (0.2e1 * (0.4e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) + g) * l1 * m2 * m2 * sin(x[0]) * d2 * d2 * cos(x[0] + x[1]) * sin(x[0] + x[1]) + l1 * m2 * m2 * (0.4e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - 0.2e1 * x[2] * l1 * (x[2] + x[3])) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + l1 * m2 * (-m2 * (0.4e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - 0.2e1 * x[2] * l1 * (x[2] + x[3])) * d2 * sin(x[0] + x[1]) + x[3] * cos(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 + x[3] * cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + (g * (d1 * m1 + l1 * m2) * sin(x[0]) + ((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3] + 0.2e1 * u[1] - u[0]) * cos(x[0]) + ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * x[2] * sin(x[0])) * d2 * sin(x[0] + x[1]) + m2 * (-x[3] * cos(x[2]) * sin(x[0]) * d2 * l1 * l1 * m2 * cos(x[2] + x[3]) - x[3] * sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * l1 * m2 + g * l1 * (d1 * m1 + l1 * m2) * pow(cos(x[0]), 0.2e1) + l1 * ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * x[2] * cos(x[0]) - 0.2e1 * l1 * (((-x[2] / 0.2e1 - x[3] / 0.2e1) * d2 * d2 + l1 * l1 * x[2] / 0.2e1) * m2 + (d1 * d1 * m1 / 0.2e1 + J1 / 0.2e1 - J2 / 0.2e1) * x[2] - J2 * x[3] / 0.2e1 + u[1] - u[0] / 0.2e1) * sin(x[0]) + g * (d1 * d1 * m1 - d1 * l1 * m1 + J1)) * d2 * cos(x[0] + x[1])) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1)) + lmd[3] * (-(0.4e1 * x[2] * l1 * (x[2] + x[3]) * cos(x[0]) + g) * l1 * m2 * m2 * sin(x[0]) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - l1 * m2 * (-m2 * (0.4e1 * x[2] * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + g * cos(x[0]) - 0.2e1 * x[2] * l1 * (x[2] + x[3])) * d2 * sin(x[0] + x[1]) + x[3] * cos(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 + x[3] * cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + (g * (d1 * m1 + l1 * m2) * sin(x[0]) + ((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3] + 0.2e1 * u[1] - u[0]) * cos(x[0]) + ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * x[2] * sin(x[0])) * d2 * cos(x[0] + x[1]) + m2 * (-x[3] * cos(x[2]) * sin(x[0]) * d2 * l1 * l1 * m2 * cos(x[2] + x[3]) - x[3] * sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * l1 * m2 + g * l1 * (d1 * m1 + l1 * m2) * pow(cos(x[0]), 0.2e1) + l1 * ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * x[2] * cos(x[0]) - 0.2e1 * l1 * (((-x[2] / 0.2e1 - x[3] / 0.2e1) * d2 * d2 + l1 * l1 * x[2] / 0.2e1) * m2 + (d1 * d1 * m1 / 0.2e1 + J1 / 0.2e1 - J2 / 0.2e1) * x[2] - J2 * x[3] / 0.2e1 + u[1] - u[0] / 0.2e1) * sin(x[0]) + g * (d1 * d1 * m1 - d1 * l1 * m1 + J1)) * d2 * sin(x[0] + x[1]) + l1 * m2 * (((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3]) * cos(x[2]) * d2 * cos(x[2] + x[3]) + l1 * sin(x[2]) * m2 * (((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3]) * d2 * sin(x[2] + x[3]) + 0.2e1 * x[2] * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * (x[2] + x[3]) * cos(x[0]) - ((d1 * d2 * d2 * m1 + J2 * l1) * m2 + J2 * d1 * m1) * g * sin(x[0]) + d2 * d2 * l1 * l1 * m2 * m2 * x[3] + (((d1 * d1 * m1 + J1) * x[3] - u[1] + u[0]) * d2 * d2 - l1 * l1 * (-J2 * x[3] + u[1])) * m2 + J2 * (d1 * d1 * m1 + J1) * x[3] + (-u[1] + u[0]) * J2 - u[1] * (d1 * d1 * m1 + J1)) * pow((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1), -0.2e1) * (-0.2e1 * (0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * cos(x[0] + x[1]) * sin(x[0] + x[1]) - 0.2e1 * cos(x[0]) * pow(sin(x[0] + x[1]), 0.2e1) * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 + 0.2e1 * cos(x[0]) * sin(x[0]) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1)) );
	hx[2] = -( -(x[2] - xf[2]) * q[2] - lmd[0] - lmd[2] * (l1 * m2 * m2 * (0.2e1 * l1 * (x[2] + x[3]) * cos(x[0]) + 0.2e1 * x[2] * l1 * cos(x[0])) * sin(x[0]) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - l1 * m2 * (m2 * (0.2e1 * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + 0.2e1 * x[2] * l1 * pow(cos(x[0]), 0.2e1) - l1 * (x[2] + x[3]) - x[2] * l1) * d2 * sin(x[0] + x[1]) + cos(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 + cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + (d2 * d2 * m2 + J2) * cos(x[0]) - sin(x[0]) * (d2 * d2 * m2 + J2) * (x[2] + x[3]) - x[2] * sin(x[0]) * (d2 * d2 * m2 + J2)) * d2 * cos(x[0] + x[1]) + l1 * m2 * (-cos(x[2]) * sin(x[0]) * d2 * l1 * m2 * cos(x[2] + x[3]) - sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 - (d2 * d2 * m2 + J2) * (x[2] + x[3]) * cos(x[0]) - x[2] * (d2 * d2 * m2 + J2) * cos(x[0]) + sin(x[0]) * (-d2 * d2 * m2 - J2)) * d2 * sin(x[0] + x[1]) + cos(x[2]) * d2 * l1 * m2 * (d2 * d2 * m2 + J2) * cos(x[2] + x[3]) + sin(x[2]) * d2 * l1 * m2 * (d2 * d2 * m2 + J2) * sin(x[2] + x[3]) - sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * (x[2] + x[3]) * cos(x[0]) - x[2] * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * cos(x[0]) + (d2 * d2 * m2 + J2) * (d1 * d1 * m1 + l1 * l1 * m2 + J1)) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1)) - lmd[3] * (-(0.4e1 * l1 * (x[2] + x[3]) * cos(x[0]) + 0.4e1 * x[2] * l1 * cos(x[0])) * l1 * m2 * m2 * sin(x[0]) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - l1 * m2 * (-m2 * (0.4e1 * l1 * (x[2] + x[3]) * pow(cos(x[0]), 0.2e1) + 0.4e1 * x[2] * l1 * pow(cos(x[0]), 0.2e1) - 0.2e1 * l1 * (x[2] + x[3]) - 0.2e1 * x[2] * l1) * d2 * sin(x[0] + x[1]) + ((-d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 - J2) * cos(x[0]) + ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * x[2] * sin(x[0]) + ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * sin(x[0])) * d2 * cos(x[0] + x[1]) + m2 * (l1 * ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * x[2] * cos(x[0]) + l1 * ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * (x[2] + x[3]) * cos(x[0]) - 0.2e1 * l1 * ((-d2 * d2 / 0.2e1 + l1 * l1 / 0.2e1) * m2 + d1 * d1 * m1 / 0.2e1 + J1 / 0.2e1 - J2 / 0.2e1) * sin(x[0])) * d2 * sin(x[0] + x[1]) + l1 * m2 * ((-d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 - J2) * cos(x[2]) * d2 * cos(x[2] + x[3]) + l1 * sin(x[2]) * m2 * ((-d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 - J2) * d2 * sin(x[2] + x[3]) + 0.2e1 * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * (x[2] + x[3]) * cos(x[0]) + 0.2e1 * x[2] * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * cos(x[0])) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1)) );
	hx[3] = -( -(x[3] - xf[3]) * q[3] - lmd[1] - lmd[2] * (0.2e1 * l1 * l1 * m2 * m2 * x[2] * cos(x[0]) * sin(x[0]) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - l1 * m2 * (m2 * (0.2e1 * x[2] * l1 * pow(cos(x[0]), 0.2e1) - x[2] * l1) * d2 * sin(x[0] + x[1]) - x[2] * cos(x[0]) * cos(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + x[2] * cos(x[0]) * sin(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 + (d2 * d2 * m2 + J2) * cos(x[0]) - x[2] * sin(x[0]) * (d2 * d2 * m2 + J2)) * d2 * cos(x[0] + x[1]) + l1 * m2 * (x[2] * cos(x[2]) * sin(x[0]) * d2 * l1 * m2 * sin(x[2] + x[3]) - x[2] * sin(x[2]) * sin(x[0]) * d2 * l1 * m2 * cos(x[2] + x[3]) - x[2] * (d2 * d2 * m2 + J2) * cos(x[0]) + sin(x[0]) * (-d2 * d2 * m2 - J2)) * d2 * sin(x[0] + x[1]) + cos(x[2]) * d2 * l1 * m2 * (d2 * d2 * m2 + J2) * cos(x[2] + x[3]) - cos(x[2]) * d2 * l1 * m2 * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * sin(x[2] + x[3]) + sin(x[2]) * d2 * l1 * m2 * (d2 * d2 * m2 + J2) * sin(x[2] + x[3]) + sin(x[2]) * d2 * l1 * m2 * (d2 * d2 * m2 + J2) * (x[2] + x[3]) * cos(x[2] + x[3]) - x[2] * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * cos(x[0])) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1)) - lmd[3] * (-0.4e1 * l1 * l1 * m2 * m2 * x[2] * cos(x[0]) * sin(x[0]) * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) - l1 * m2 * (-m2 * (0.4e1 * x[2] * l1 * pow(cos(x[0]), 0.2e1) - 0.2e1 * x[2] * l1) * d2 * sin(x[0] + x[1]) + cos(x[0]) * cos(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 - x[3] * cos(x[0]) * cos(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + cos(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * m2 + x[3] * cos(x[0]) * sin(x[2]) * cos(x[2] + x[3]) * d2 * l1 * m2 + (-d2 * d2 * m2 - J2) * cos(x[0]) + ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * x[2] * sin(x[0])) * d2 * cos(x[0] + x[1]) + m2 * (-cos(x[2]) * sin(x[0]) * d2 * l1 * l1 * m2 * cos(x[2] + x[3]) + x[3] * cos(x[2]) * sin(x[0]) * d2 * l1 * l1 * m2 * sin(x[2] + x[3]) - sin(x[0]) * sin(x[2]) * sin(x[2] + x[3]) * d2 * l1 * l1 * m2 - x[3] * sin(x[2]) * sin(x[0]) * d2 * l1 * l1 * m2 * cos(x[2] + x[3]) + l1 * ((d2 * d2 + l1 * l1) * m2 + d1 * d1 * m1 + J1 + J2) * x[2] * cos(x[0]) - 0.2e1 * l1 * (-d2 * d2 * m2 / 0.2e1 - J2 / 0.2e1) * sin(x[0])) * d2 * sin(x[0] + x[1]) + l1 * m2 * (-d2 * d2 * m2 - J2) * cos(x[2]) * d2 * cos(x[2] + x[3]) - l1 * m2 * (((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3]) * cos(x[2]) * d2 * sin(x[2] + x[3]) + l1 * sin(x[2]) * m2 * (-d2 * d2 * m2 - J2) * d2 * sin(x[2] + x[3]) + l1 * m2 * (((-x[2] - x[3]) * d2 * d2 + l1 * l1 * x[2]) * m2 + (d1 * d1 * m1 + J1 - J2) * x[2] - J2 * x[3]) * sin(x[2]) * d2 * cos(x[2] + x[3]) + 0.2e1 * x[2] * sin(x[0]) * d2 * d2 * l1 * l1 * m2 * m2 * cos(x[0]) + d2 * d2 * l1 * l1 * m2 * m2 + ((d1 * d1 * m1 + J1) * d2 * d2 + J2 * l1 * l1) * m2 + J2 * (d1 * d1 * m1 + J1)) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1)) );
}

void nmpc_model::hufunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd& lmd, Eigen::Ref<Eigen::VectorXd> hu)
{
	hu[0] = u[0] * r[0] + lmd[2] * (-d2 * d2 * m2 - J2) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1)) + lmd[3] * (l1 * m2 * cos(x[0]) * d2 * cos(x[0] + x[1]) + m2 * l1 * sin(x[0]) * d2 * sin(x[0] + x[1]) + d2 * d2 * m2 + J2) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1));
	hu[1] = u[1] * r[1] + lmd[2] * (l1 * m2 * cos(x[0]) * d2 * cos(x[0] + x[1]) + m2 * l1 * sin(x[0]) * d2 * sin(x[0] + x[1]) + d2 * d2 * m2 + J2) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1)) + lmd[3] * (-0.2e1 * l1 * m2 * cos(x[0]) * d2 * cos(x[0] + x[1]) - 0.2e1 * m2 * l1 * sin(x[0]) * d2 * sin(x[0] + x[1]) + (-d2 * d2 - l1 * l1) * m2 - J2 - d1 * d1 * m1 - J1) / ((0.2e1 * pow(cos(x[0]), 0.2e1) - 0.1e1) * l1 * l1 * m2 * m2 * d2 * d2 * pow(cos(x[0] + x[1]), 0.2e1) + 0.2e1 * cos(x[0]) * cos(x[0] + x[1]) * sin(x[0]) * sin(x[0] + x[1]) * d2 * d2 * l1 * l1 * m2 * m2 - pow(cos(x[0]), 0.2e1) * d2 * d2 * l1 * l1 * m2 * m2 + ((-d1 * d1 * m1 - J1) * d2 * d2 - J2 * l1 * l1) * m2 - J2 * (d1 * d1 * m1 + J1));
}


double nmpc_model::stagecost(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u)
{
	return 0.5 * (q[0] * (x[0] - xf[0]) * (x[0] - xf[0]) + q[1] * (x[1] - xf[1]) * (x[1] - xf[1]) + q[2] * (x[2] - xf[2]) * (x[2] - xf[2]) + q[3] * (x[3] - xf[3]) * (x[3] - xf[3]) + r[0] * u[0] * u[0] + r[1] * u[1] * u[1]);
}

double nmpc_model::terminalcost(const double t, const Eigen::VectorXd& x)
{
	return 0.5 * (sf[0] * (x[0] - xf[0]) * (x[0] - xf[0]) + sf[1] * (x[1] - xf[1]) * (x[1] - xf[1]) + sf[2] * (x[2] - xf[2]) * (x[2] - xf[2]) + sf[3] * (x[3] - xf[3]) * (x[3] - xf[3]));
}

void nmpc_model::Lufunc(const double t, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> Ju)
{
	Ju[0] = r[0] * u[0];
	Ju[1] = r[1] * u[1];
}