#include"fvm_2d.h"
#include"output.h"
#include<cmath>
void output2d_blow_up(Fluid2d* fluids, Block2d block);
int gausspoint = 0; //initilization
double* gauss_loc = new double[gausspoint]; //initilization
double* gauss_weight = new double[gausspoint]; //initilization
GKS2d_type gks2dsolver = kfvs1st_2d; //initilization

Reconstruction_within_Cell_2D_normal cellreconstruction_2D_normal = WENO5_AO_normal; //initilization
Reconstruction_within_Cell_2D_normal_du cellreconstruction_2D_normal_du = WENO5_normal_du; //initilization

Reconstruction_within_Cell_2D_tangent cellreconstruction_2D_tangent = WENO5_AO_tangent; //initilization
Reconstruction_forG0_2D_normal g0reconstruction_2D_normal = Center_do_nothing_normal; //initilization
Reconstruction_forG0_2D_tangent g0reconstruction_2D_tangent = Center_all_collision_multi; //initilization

Flux_function_2d flux_function_2d = CGKS2D; //initilization
TimeMarchingCoefficient_2d timecoe_list_2d = S1O1_2D; //initilization




void SetGuassPoint()
{
	if (gausspoint == 0)
	{
		cout << "no gausspoint specify" << endl;
		exit(0);
	}
	if (gausspoint == 1)
	{
		gauss_loc = new double[1]; gauss_loc[0] = 0.0;
		gauss_weight = new double[1]; gauss_weight[0] = 1.0;
	}
	if (gausspoint == 2)
	{
		gauss_loc = new double[2]; gauss_loc[0] = -sqrt(1.0 / 3.0); gauss_loc[1] = -gauss_loc[0];
		gauss_weight = new double[2]; gauss_weight[0] = 0.5; gauss_weight[1] = 0.5;
	}
	if (gausspoint == 3)
	{
		cout << "the usage of 3 gausspoint may be imcompatiable for some reconstructions" << endl;
		//exit(0);
		gauss_loc = new double[3]; gauss_loc[0] = 0.0; gauss_loc[1] = sqrt(3.0 / 5.0); gauss_loc[2] = -gauss_loc[1];
		gauss_weight = new double[3];
		gauss_weight[0] = 4.0 / 9.0; gauss_weight[1] = 0.5 - 0.5 * gauss_weight[0]; gauss_weight[2] = gauss_weight[1];
	}
	if (gausspoint == 4)
	{
		cout << "the usage of 4 gausspoint may be imcompatiable for some reconstructions" << endl;
		gauss_loc = new double[4]; gauss_weight = new double[4];

		gauss_loc[0] = -1.0; gauss_loc[1] = -gauss_loc[0];
		gauss_loc[2] = -sqrt(5) / 5.0; gauss_loc[3] = -gauss_loc[2];
		gauss_weight[0] = 1.0 / 12.0;		gauss_weight[1] = 1.0 / 12.0;
		gauss_weight[2] = 5.0 / 12.0;		gauss_weight[3] = 5.0 / 12.0;

	}
	cout << gausspoint << " gausspoint(s) specify" << endl;
}


double Get_CFL1(Block2d& block, Fluid2d* fluids, double tstop)
{
	double dt = 111.0, dtmin = 111.0;

	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			dt = Dtx(dt, fluids[i * (block.ny) + j].dx, block.CFL, fluids[i * (block.ny) + j].primvar[0], fluids[i * (block.ny) + j].primvar[1],
				fluids[i * (block.ny) + j].primvar[2], fluids[i * (block.ny) + j].primvar[3]);

			if (dt < dtmin)
			{
				dtmin = dt;
			}

		}
	}
	dt = dtmin;

	if (block.step % 100 == 0)
	{
		cout << "step= " << block.step
			<< "time size is " << dt
			<< " time= " << block.t
			<< endl;
	}
	return dt;
}

double Thetamax(Interface2d* xinterfaces, Interface2d* yinterfaces, Block2d& block, Fluid2d* fluids)
{
	double thetamax = 0.0;

	for (int i = block.ghost - 1; i < block.nx - block.ghost + 1; i++)
	{
		for (int j = block.ghost - 1; j < block.ny - block.ghost + 1; j++)
		{
			int index = i * block.ny + j;
			double thetamax_x = max(xinterfaces[i * (block.ny + 1) + j].line.right.gama, xinterfaces[i * (block.ny + 1) + j].line.left.gama);
			double thetamax_y = max(yinterfaces[i * (block.ny + 1) + j].line.right.gama, yinterfaces[i * (block.ny + 1) + j].line.left.gama);

		

			thetamax = max({thetamax, thetamax_x, thetamax_y});
		}
	}

	//cout << block.thetamax << endl;

	return thetamax;
}
void A_point(double* a, double* der, double* prim)
{
	double R4, R3, R2;
	double overden = 1.0 / prim[0];
	R4 = der[3] * overden - 0.5 * (prim[1] * prim[1] + prim[2] * prim[2] + 0.5 * (K + 2) / prim[3]) * der[0] * overden;
	R3 = (der[2] - prim[2] * der[0]) * overden;
	R2 = (der[1] - prim[1] * der[0]) * overden;
	a[3] = (4.0 / (K + 2)) * prim[3] * prim[3] * (2 * R4 - 2 * prim[1] * R2 - 2 * prim[2] * R3);
	a[2] = 2 * prim[3] * R3 - prim[2] * a[3];
	a[1] = 2 * prim[3] * R2 - prim[1] * a[3];
	a[0] = der[0] * overden - prim[1] * a[1] - prim[2] * a[2] - 0.5 * a[3] * (prim[1] * prim[1] + prim[2] * prim[2] + 0.5 * (K + 2) / prim[3]);
}
void G_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF& m)
{

	psi[0] = a[0] * m.uvxi[no_u][no_v][no_xi] + a[1] * m.uvxi[no_u + 1][no_v][no_xi] + a[2] * m.uvxi[no_u][no_v + 1][no_xi] + a[3] * 0.5 * (m.uvxi[no_u + 2][no_v][no_xi] + m.uvxi[no_u][no_v + 2][no_xi] + m.uvxi[no_u][no_v][no_xi + 1]);
	psi[1] = a[0] * m.uvxi[no_u + 1][no_v][no_xi] + a[1] * m.uvxi[no_u + 2][no_v][no_xi] + a[2] * m.uvxi[no_u + 1][no_v + 1][no_xi] + a[3] * 0.5 * (m.uvxi[no_u + 3][no_v][no_xi] + m.uvxi[no_u + 1][no_v + 2][no_xi] + m.uvxi[no_u + 1][no_v][no_xi + 1]);
	psi[2] = a[0] * m.uvxi[no_u][no_v + 1][no_xi] + a[1] * m.uvxi[no_u + 1][no_v + 1][no_xi] + a[2] * m.uvxi[no_u][no_v + 2][no_xi] + a[3] * 0.5 * (m.uvxi[no_u + 2][no_v + 1][no_xi] + m.uvxi[no_u][no_v + 3][no_xi] + m.uvxi[no_u][no_v + 1][no_xi + 1]);
	psi[3] = 0.5 * (a[0] * (m.uvxi[no_u + 2][no_v][no_xi] + m.uvxi[no_u][no_v + 2][no_xi] + m.uvxi[no_u][no_v][no_xi + 1]) +
		a[1] * (m.uvxi[no_u + 3][no_v][no_xi] + m.uvxi[no_u + 1][no_v + 2][no_xi] + m.uvxi[no_u + 1][no_v][no_xi + 1]) +
		a[2] * (m.uvxi[no_u + 2][no_v + 1][no_xi] + m.uvxi[no_u][no_v + 3][no_xi] + m.uvxi[no_u][no_v + 1][no_xi + 1]) +
		a[3] * 0.5 * (m.uvxi[no_u + 4][no_v][no_xi] + m.uvxi[no_u][no_v + 4][no_xi] + m.uvxi[no_u][no_v][no_xi + 2] + 2 * m.uvxi[no_u + 2][no_v + 2][no_xi] + 2 * m.uvxi[no_u + 2][no_v][no_xi + 1] + 2 * m.uvxi[no_u][no_v + 2][no_xi + 1]));
}
void GL_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF& m)
{

	psi[0] = a[0] * m.upvxi[no_u][no_v][no_xi] + a[1] * m.upvxi[no_u + 1][no_v][no_xi] + a[2] * m.upvxi[no_u][no_v + 1][no_xi] + a[3] * 0.5 * (m.upvxi[no_u + 2][no_v][no_xi] + m.upvxi[no_u][no_v + 2][no_xi] + m.upvxi[no_u][no_v][no_xi + 1]);
	psi[1] = a[0] * m.upvxi[no_u + 1][no_v][no_xi] + a[1] * m.upvxi[no_u + 2][no_v][no_xi] + a[2] * m.upvxi[no_u + 1][no_v + 1][no_xi] + a[3] * 0.5 * (m.upvxi[no_u + 3][no_v][no_xi] + m.upvxi[no_u + 1][no_v + 2][no_xi] + m.upvxi[no_u + 1][no_v][no_xi + 1]);
	psi[2] = a[0] * m.upvxi[no_u][no_v + 1][no_xi] + a[1] * m.upvxi[no_u + 1][no_v + 1][no_xi] + a[2] * m.upvxi[no_u][no_v + 2][no_xi] + a[3] * 0.5 * (m.upvxi[no_u + 2][no_v + 1][no_xi] + m.upvxi[no_u][no_v + 3][no_xi] + m.upvxi[no_u][no_v + 1][no_xi + 1]);
	psi[3] = 0.5 * (a[0] * (m.upvxi[no_u + 2][no_v][no_xi] + m.upvxi[no_u][no_v + 2][no_xi] + m.upvxi[no_u][no_v][no_xi + 1]) +
		a[1] * (m.upvxi[no_u + 3][no_v][no_xi] + m.upvxi[no_u + 1][no_v + 2][no_xi] + m.upvxi[no_u + 1][no_v][no_xi + 1]) +
		a[2] * (m.upvxi[no_u + 2][no_v + 1][no_xi] + m.upvxi[no_u][no_v + 3][no_xi] + m.upvxi[no_u][no_v + 1][no_xi + 1]) +
		a[3] * 0.5 * (m.upvxi[no_u + 4][no_v][no_xi] + m.upvxi[no_u][no_v + 4][no_xi] + m.upvxi[no_u][no_v][no_xi + 2] + 2 * m.upvxi[no_u + 2][no_v + 2][no_xi] + 2 * m.upvxi[no_u + 2][no_v][no_xi + 1] + 2 * m.upvxi[no_u][no_v + 2][no_xi + 1]));
}
void GR_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF& m)
{
	// Similar to the GL_address
	psi[0] = a[0] * m.unvxi[no_u][no_v][no_xi] + a[1] * m.unvxi[no_u + 1][no_v][no_xi] + a[2] * m.unvxi[no_u][no_v + 1][no_xi] + a[3] * 0.5 * (m.unvxi[no_u + 2][no_v][no_xi] + m.unvxi[no_u][no_v + 2][no_xi] + m.unvxi[no_u][no_v][no_xi + 1]);
	psi[1] = a[0] * m.unvxi[no_u + 1][no_v][no_xi] + a[1] * m.unvxi[no_u + 2][no_v][no_xi] + a[2] * m.unvxi[no_u + 1][no_v + 1][no_xi] + a[3] * 0.5 * (m.unvxi[no_u + 3][no_v][no_xi] + m.unvxi[no_u + 1][no_v + 2][no_xi] + m.unvxi[no_u + 1][no_v][no_xi + 1]);
	psi[2] = a[0] * m.unvxi[no_u][no_v + 1][no_xi] + a[1] * m.unvxi[no_u + 1][no_v + 1][no_xi] + a[2] * m.unvxi[no_u][no_v + 2][no_xi] + a[3] * 0.5 * (m.unvxi[no_u + 2][no_v + 1][no_xi] + m.unvxi[no_u][no_v + 3][no_xi] + m.unvxi[no_u][no_v + 1][no_xi + 1]);
	psi[3] = 0.5 * (a[0] * (m.unvxi[no_u + 2][no_v][no_xi] + m.unvxi[no_u][no_v + 2][no_xi] + m.unvxi[no_u][no_v][no_xi + 1]) +
		a[1] * (m.unvxi[no_u + 3][no_v][no_xi] + m.unvxi[no_u + 1][no_v + 2][no_xi] + m.unvxi[no_u + 1][no_v][no_xi + 1]) +
		a[2] * (m.unvxi[no_u + 2][no_v + 1][no_xi] + m.unvxi[no_u][no_v + 3][no_xi] + m.unvxi[no_u][no_v + 1][no_xi + 1]) +
		a[3] * 0.5 * (m.unvxi[no_u + 4][no_v][no_xi] + m.unvxi[no_u][no_v + 4][no_xi] + m.unvxi[no_u][no_v][no_xi + 2] + 2 * m.unvxi[no_u + 2][no_v + 2][no_xi] + 2 * m.unvxi[no_u + 2][no_v][no_xi + 1] + 2 * m.unvxi[no_u][no_v + 2][no_xi + 1]));
}


void ICfor_uniform_2d(Fluid2d* fluid, double* prim, Block2d block)
{

	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			fluid[i * block.ny + j].primvar[0] = prim[0];
			fluid[i * block.ny + j].primvar[1] = prim[1];
			fluid[i * block.ny + j].primvar[2] = prim[2];
			fluid[i * block.ny + j].primvar[3] = prim[3];
		}

	}

	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			Primvar_to_convar_2D(fluid[i * block.ny + j].convar, fluid[i * block.ny + j].primvar);

		}
	}
}

bool negative_density_or_pressure(double* primvar)
{
	//detect whether density or pressure is negative
	if (primvar[0] < 0 || primvar[3] < 0 ||
		isnan(primvar[0])
		|| isnan(primvar[3]))
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Convar_to_Primvar(Fluid2d* fluids, Block2d block)
{
	bool all_positive = true;
#pragma omp parallel  for
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			Convar_to_primvar_2D(fluids[i * block.ny + j].primvar, fluids[i * block.ny + j].convar);
			if (negative_density_or_pressure(fluids[i * block.ny + j].primvar))
			{
				all_positive = false;
			}
		}
	}

	if (all_positive == false)
	{
		cout << "the program blows up at t=" << block.t << "!" << endl;
		output2d_blow_up(fluids, block);
		exit(0);
	}
}
void Convar_to_primvar(Fluid2d* fluid, Block2d block)
{
#pragma omp parallel  for
	for (int i = 0; i < block.nx * block.ny; i++)
	{
		Convar_to_primvar_2D(fluid[i].primvar, fluid[i].convar);
	}
}
void YchangetoX(double* fluidtmp, double* fluid)
{
	fluidtmp[0] = fluid[0];
	fluidtmp[1] = fluid[2];
	fluidtmp[2] = -fluid[1];
	fluidtmp[3] = fluid[3];
}
void XchangetoY(double* fluidtmp, double* fluid)
{
	fluidtmp[0] = fluid[0];
	fluidtmp[1] = -fluid[2];
	fluidtmp[2] = fluid[1];
	fluidtmp[3] = fluid[3];
}
void CopyFluid_new_to_old(Fluid2d* fluids, Block2d block)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.ghost + block.nodex; i++)
	{
		for (int j = block.ghost; j < block.ghost + block.nodey; j++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].convar_old[var] = fluids[i * block.ny + j].convar[var];
			}
		}
	}
}

void Reconstruction_within_cell(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
#pragma omp parallel for collapse(2)
	for (int i = block.ghost - 2; i < block.nx - block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.ny - block.ghost + 2; j++)
		{
			cellreconstruction_2D_normal
			(xinterfaces[i * (block.ny + 1) + j], xinterfaces[(i + 1) * (block.ny + 1) + j],
				yinterfaces[i * (block.ny + 1) + j], yinterfaces[i * (block.ny + 1) + j + 1], &fluids[i * block.ny + j], block);
		}
	}

#pragma omp parallel for collapse(2)
	for (int i = block.ghost - 2; i < block.nx - block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.ny - block.ghost + 2; j++)
		{
			for (int var = 0; var < 4; var++)
			{
				xinterfaces[(i + 1) * (block.ny + 1) + j].gauss[0].right.convar[var] = xinterfaces[(i + 1) * (block.ny + 1) + j].line.right.convar[var];
				xinterfaces[(i + 1) * (block.ny + 1) + j].gauss[0].left.convar[var] = xinterfaces[(i + 1) * (block.ny + 1) + j].line.left.convar[var];
				yinterfaces[i * (block.ny + 1) + j].gauss[0].right.convar[var] = yinterfaces[i * (block.ny + 1) + j].line.right.convar[var];
				yinterfaces[i * (block.ny + 1) + j + 1].gauss[0].left.convar[var] = yinterfaces[i * (block.ny + 1) + j + 1].line.left.convar[var];
			}
		}
	}
}

void Reconstruction_within_cell0(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block, int stage)
{
#pragma omp parallel for collapse(2)
	for (int i = block.ghost - 2; i < block.nx - block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.ny - block.ghost + 2; j++)
		{
			TENO70_normal
			(xinterfaces[i * (block.ny + 1) + j], xinterfaces[(i + 1) * (block.ny + 1) + j],
				yinterfaces[i * (block.ny + 1) + j], yinterfaces[i * (block.ny + 1) + j + 1], &fluids[i * block.ny + j], block);
		}
	}

#pragma omp parallel for collapse(2)
	for (int i = block.ghost - 2; i < block.nx - block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.ny - block.ghost + 2; j++)
		{
			for (int var = 0; var < 4; var++)
			{
				xinterfaces[(i + 1) * (block.ny + 1) + j].gauss[0].right.convar[var] = xinterfaces[(i + 1) * (block.ny + 1) + j].line.right.convar[var];
				xinterfaces[(i + 1) * (block.ny + 1) + j].gauss[0].left.convar[var] = xinterfaces[(i + 1) * (block.ny + 1) + j].line.left.convar[var];
				yinterfaces[i * (block.ny + 1) + j].gauss[0].right.convar[var] = yinterfaces[i * (block.ny + 1) + j].line.right.convar[var];
				yinterfaces[i * (block.ny + 1) + j + 1].gauss[0].left.convar[var] = yinterfaces[i * (block.ny + 1) + j + 1].line.left.convar[var];
			}

			xinterfaces[(i + 1) * (block.ny + 1) + j].line.right.gama = xinterfaces[(i + 1) * (block.ny + 1) + j].line.right.der1x[0];
			xinterfaces[(i + 1) * (block.ny + 1) + j].line.left.gama = xinterfaces[(i + 1) * (block.ny + 1) + j].line.left.der1x[0];
			yinterfaces[i * (block.ny + 1) + j + 1].line.right.gama = yinterfaces[i * (block.ny + 1) + j + 1].line.right.der1x[0];
			yinterfaces[i * (block.ny + 1) + j + 1].line.left.gama = yinterfaces[i * (block.ny + 1) + j + 1].line.left.der1x[0];

		}
	}

}

void Reconstruction_within_cell1(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
#pragma omp parallel for collapse(2)
	for (int i = block.ghost - 2; i < block.nx - block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.ny - block.ghost + 2; j++)
		{
				cellreconstruction_2D_normal
				(xinterfaces[i * (block.ny + 1) + j], xinterfaces[(i + 1) * (block.ny + 1) + j],
					yinterfaces[i * (block.ny + 1) + j], yinterfaces[i * (block.ny + 1) + j + 1], &fluids[i * block.ny + j], block);
		}
	}

#pragma omp parallel for collapse(2)
	for (int i = block.ghost - 2; i < block.nx - block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.ny - block.ghost + 2; j++)
		{
			int index = i * block.ny + j;

			for (int var = 0; var < 4; var++)
			{
				xinterfaces[(i + 1) * (block.ny + 1) + j].gauss[0].right.convar[var] = xinterfaces[(i + 1) * (block.ny + 1) + j].line.right.convar[var];
				xinterfaces[(i + 1) * (block.ny + 1) + j].gauss[0].left.convar[var] = xinterfaces[(i + 1) * (block.ny + 1) + j].line.left.convar[var];
				yinterfaces[i * (block.ny + 1) + j].gauss[0].right.convar[var] = yinterfaces[i * (block.ny + 1) + j].line.right.convar[var];
				yinterfaces[i * (block.ny + 1) + j + 1].gauss[0].left.convar[var] = yinterfaces[i * (block.ny + 1) + j + 1].line.left.convar[var];
			}

	
			fluids[index].sensor = max({ xinterfaces[i * (block.ny + 1) + j].line.right.sensor, xinterfaces[i * (block.ny + 1) + j].line.left.sensor,
				                         yinterfaces[i * (block.ny + 1) + j].line.right.sensor, yinterfaces[i * (block.ny + 1) + j].line.left.sensor }); //cell max
		}
	}


}


void TENO7origin_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)
{

	if (block.uniform == true)
	{
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 2))
		{
			TENO7origin(left.line.right, right.line.left,
				fluids[-3 * block.ny].convar, fluids[-2 * block.ny].convar, fluids[-block.ny].convar,
				fluids[0].convar,
				fluids[block.ny].convar, fluids[2 * block.ny].convar, fluids[3 * block.ny].convar,
				fluids[0].dx);
		}

		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 2))
		{
			double wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4];
			YchangetoX(wn3tmp, fluids[-3].convar);
			YchangetoX(wn2tmp, fluids[-2].convar);
			YchangetoX(wn1tmp, fluids[-1].convar);
			YchangetoX(wtmp, fluids[0].convar);
			YchangetoX(wp1tmp, fluids[1].convar);
			YchangetoX(wp2tmp, fluids[2].convar);
			YchangetoX(wp3tmp, fluids[3].convar);

			TENO7origin(down.line.right, up.line.left,
				wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp,
				fluids[0].dy);
		}
	}
	else
	{
		// for non-uniform mesh.
		Point2d voidpoint;
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 2))
		{
			double dx = fluids[0].dx;
			double normal[2];
			double wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4];

			// cell left reconstruction
			Copy_Array(normal, left.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3 * block.ny].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2 * block.ny].convar, normal);
			Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);
			Global_to_Local(wp2tmp, fluids[2 * block.ny].convar, normal);
			Global_to_Local(wp3tmp, fluids[3 * block.ny].convar, normal);
			TENO7origin(left.line.right, voidpoint, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dx);

			// cell right reconstruction
			Copy_Array(normal, right.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3 * block.ny].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2 * block.ny].convar, normal);
			Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);
			Global_to_Local(wp2tmp, fluids[2 * block.ny].convar, normal);
			Global_to_Local(wp3tmp, fluids[3 * block.ny].convar, normal);
			TENO7origin(voidpoint, right.line.left, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dx);
		}

		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 2))
		{
			double dy = fluids[0].dy;
			double normal[2];
			double wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4];

			// down interface reconstruction
			Copy_Array(normal, down.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2].convar, normal);
			Global_to_Local(wn1tmp, fluids[-1].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[1].convar, normal);
			Global_to_Local(wp2tmp, fluids[2].convar, normal);
			Global_to_Local(wp3tmp, fluids[3].convar, normal);
			TENO7origin(down.line.right, voidpoint, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dy);

			// up interface reconstruction
			Copy_Array(normal, up.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2].convar, normal);
			Global_to_Local(wn1tmp, fluids[-1].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[1].convar, normal);
			Global_to_Local(wp2tmp, fluids[2].convar, normal);
			Global_to_Local(wp3tmp, fluids[3].convar, normal);
			TENO7origin(voidpoint, up.line.left, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dy);
		}
	}
}
void TENO7origin(Point2d& left, Point2d& right, double* wn3, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double* wp3, double h)
{
	double ren3[4], ren2[4], ren1[4], re0[4], rep1[4], rep2[4], rep3[4];
	double var[4], der1[4], der2[4];

	double base_left[4];
	double base_right[4];
	double wn3_primvar[4], wn2_primvar[4], wn1_primvar[4], w_primvar[4], wp1_primvar[4], wp2_primvar[4], wp3_primvar[4];
	//Convar_to_primvar_2D(wn3_primvar, wn3);
	//Convar_to_primvar_2D(wn2_primvar, wn2);
	Convar_to_primvar_2D(wn1_primvar, wn1);
	Convar_to_primvar_2D(w_primvar, w);
	Convar_to_primvar_2D(wp1_primvar, wp1);
	//Convar_to_primvar_2D(wp2_primvar, wp2);
	//Convar_to_primvar_2D(wp3_primvar, wp3);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
		base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren3[i] = wn3[i];
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
			rep3[i] = wp3[i];
		}
	}
	else
	{
		Convar_to_char(ren3, base_left, wn3);
		Convar_to_char(ren2, base_left, wn2);
		Convar_to_char(ren1, base_left, wn1);
		Convar_to_char(re0, base_left, w);
		Convar_to_char(rep1, base_left, wp1);
		Convar_to_char(rep2, base_left, wp2);
		Convar_to_char(rep3, base_left, wp3);
	}


	for (int i = 0; i < 4; i++)
	{
		TENO7origin_left(var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			left.convar[i] = var[i];
			left.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(left.convar, base_left, var);
		Char_to_convar(left.der1x, base_left, der1);
	}



	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren3[i] = wn3[i];
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
			rep3[i] = wp3[i];
		}
	}
	else
	{
		Convar_to_char(ren3, base_right, wn3);
		Convar_to_char(ren2, base_right, wn2);
		Convar_to_char(ren1, base_right, wn1);
		Convar_to_char(re0, base_right, w);
		Convar_to_char(rep1, base_right, wp1);
		Convar_to_char(rep2, base_right, wp2);
		Convar_to_char(rep3, base_right, wp3);
	}


	for (int i = 0; i < 4; i++)
	{
		TENO7origin_right(var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			right.convar[i] = var[i];
			right.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(right.convar, base_right, var);
		Char_to_convar(right.der1x, base_right, der1);
	}

	Check_Order_Reduce_by_Lambda_2D(right.is_reduce_order, right.convar);
	Check_Order_Reduce_by_Lambda_2D(left.is_reduce_order, left.convar);

	if (left.is_reduce_order == true || right.is_reduce_order == true)
	{
		if (is_reduce_order_warning == true)
			cout << " TENO6-cell-splitting order reduce" << endl;
		for (int m = 0; m < 4; m++)
		{
			right.convar[m] = w[m];
			left.convar[m] = w[m];
			right.der1x[m] = 0.0;
			left.der1x[m] = 0.0;
		}
	}
}
void TENO7origin_left(double& var, double& der1, double& der2, double wn3, double wn2, double wn1, double w, double wp1, double wp2, double wp3, double h)
{
	double qleft[5];
	double dleft[5] = { 18.0 / 35.0, 3.0 / 35.0, 9.0 / 35.0, 1.0 / 35.0 , 4.0 / 35.0 };

	qleft[0] = (2.0 / 6.0) * wn1 + (5.0 / 6.0) * w + (-1.0 / 6.0) * wp1;
	qleft[1] = (11.0 / 6.0) * w + (-7.0 / 6.0) * wp1 + (2.0 / 6.0) * wp2;
	qleft[2] = (-1.0 / 6.0) * wn2 + (5.0 / 6.0) * wn1 + (2.0 / 6.0) * w;
	qleft[3] = (25.0 / 12.0) * w + (-23.0 / 12.0) * wp1 + (13.0 / 12.0) * wp2 + (-3.0 / 12.0) * wp3;
	qleft[4] = (1.0 / 12.0) * wn3 + (-5.0 / 12.0) * wn2 + (13.0 / 12.0) * wn1 + (3.0 / 12.0) * w;


	double beta[5];

	beta[0] = (13.0 / 12.0) * pow(wn1 - 2.0 * w + wp1, 2)
		+ 0.25 * pow(wn1 - wp1, 2);

	beta[1] = (13.0 / 12.0) * pow(w - 2.0 * wp1 + wp2, 2)
		+ 0.25 * pow(3.0 * w - 4.0 * wp1 + wp2, 2);

	beta[2] = (13.0 / 12.0) * pow(wn2 - 2.0 * wn1 + w, 2)
		+ 0.25 * pow(wn2 - 4.0 * wn1 + 3.0 * w, 2);

	beta[3] = (2107.0 * w * w - 9402.0 * w * wp1 + 7042.0 * w * wp2 - 1854.0 * w * wp3
		+ 11003.0 * wp1 * wp1 - 17246.0 * wp1 * wp2 + 4642.0 * wp1 * wp3
		+ 7043.0 * wp2 * wp2 - 3882.0 * wp2 * wp3
		+ 547.0 * wp3 * wp3) / 240.0;

	beta[4] = (547.0 * wn3 * wn3 - 3882.0 * wn3 * wn2 + 4642.0 * wn3 * wn1 - 1854.0 * wn3 * w
		+ 7043.0 * wn2 * wn2 - 17246.0 * wn2 * wn1 + 7042.0 * wn2 * w
		+ 11003.0 * wn1 * wn1 - 9402.0 * wn1 * w
		+ 2107.0 * w * w) / 240.0;


	double epsilon2 = 1e-40;

	double alphaleft[5];

	for (int i = 0; i < 5; i++)
	{
		if (beta[i] < 1e-30) beta[i] = 0.0;
	}


	double gamaleft[5];
	double tau7 = fabs(beta[4] - beta[3]);

	for (int i = 0; i < 5; i++)
	{
		gamaleft[i] = pow((1.0 + tau7 / (beta[i] + epsilon2)), 6);
	}
	double gamal = 0.0;
	for (int i = 0; i < 5; i++)
	{
		gamal += gamaleft[i];
	}
	double faileft[5];
	for (int i = 0; i < 5; i++)
	{
		faileft[i] = gamaleft[i] / gamal;
	}
	for (int i = 0; i < 5; i++)
	{
		if (faileft[i] < 1e-7)
		{
			alphaleft[i] = 0.0;
		}
		else
		{
			alphaleft[i] = dleft[i];
		}
	}


	double alphal = 0.0;
	for (int i = 0; i < 5; i++) alphal += alphaleft[i];

	double omegaleft[5];
	for (int i = 0; i < 5; i++) omegaleft[i] = alphaleft[i] / alphal;

	double left = 0.0;
	for (int i = 0; i < 5; i++) left += omegaleft[i] * qleft[i];


	var = left;
	der2 = 0.0;
	der1 = 0.0;
}
void TENO7origin_right(double& var, double& der1, double& der2, double wn3, double wn2, double wn1, double w, double wp1, double wp2, double wp3, double h)
{
	double qright[5];
	double dright[5] = { 18.0 / 35.0, 9.0 / 35.0, 3.0 / 35.0, 4.0 / 35.0 , 1.0 / 35.0 };

	qright[0] = (-1.0 / 6.0) * wn1 + (5.0 / 6.0) * w + (2.0 / 6.0) * wp1;
	qright[1] = (2.0 / 6.0) * w + (5.0 / 6.0) * wp1 - (1.0 / 6.0) * wp2;
	qright[2] = (2.0 / 6.0) * wn2 - (7.0 / 6.0) * wn1 + (11.0 / 6.0) * w;
	qright[3] = (3.0 / 12.0) * w + (13.0 / 12.0) * wp1 - (5.0 / 12.0) * wp2 + (1.0 / 12.0) * wp3;
	qright[4] = (-3.0 / 12.0) * wn3 + (13.0 / 12.0) * wn2 - (23.0 / 12.0) * wn1 + (25.0 / 12.0) * w;


	double beta[5];

	beta[0] = (13.0 / 12.0) * pow(wn1 - 2.0 * w + wp1, 2)
		+ 0.25 * pow(wn1 - wp1, 2);

	beta[1] = (13.0 / 12.0) * pow(w - 2.0 * wp1 + wp2, 2)
		+ 0.25 * pow(3.0 * w - 4.0 * wp1 + wp2, 2);

	beta[2] = (13.0 / 12.0) * pow(wn2 - 2.0 * wn1 + w, 2)
		+ 0.25 * pow(wn2 - 4.0 * wn1 + 3.0 * w, 2);


	beta[3] = (1.0 / 36.0) * pow(-11.0 * w + 18.0 * wp1 - 9.0 * wp2 + 2 * wp3, 2) + (13.0 / 12.0) * pow(2.0 * w - 5.0 * wp1 + 4.0 * wp2 - wp3, 2)
		+ (781.0 / 720.0) * pow(-w + 3.0 * wp1 - 3.0 * wp2 + wp3, 2);

	beta[4] = (547.0 * wn3 * wn3 - 3882.0 * wn3 * wn2 + 4642.0 * wn3 * wn1 - 1854.0 * wn3 * w
		+ 7043.0 * wn2 * wn2 - 17246.0 * wn2 * wn1 + 7042.0 * wn2 * w
		+ 11003.0 * wn1 * wn1 - 9402.0 * wn1 * w
		+ 2107.0 * w * w) / 240.0;



	for (int i = 0; i < 5; i++)
	{
		if (beta[i] < 1e-30) beta[i] = 0.0;
	}

	double epsilon2 = 1e-40;


	double alpharight[5];
	double gamaright[5];

	double tau7 = fabs(beta[4] - beta[3]);

	for (int i = 0; i < 5; i++)
	{
		gamaright[i] = pow((1.0 + tau7 / (beta[i] + epsilon2)), 6);
	}
	double gamar = 0.0;
	for (int i = 0; i < 5; i++)
	{
		gamar += gamaright[i];
	}
	double fairight[5];
	for (int i = 0; i < 5; i++)
	{
		fairight[i] = gamaright[i] / gamar;
	}
	for (int i = 0; i < 5; i++)
	{
		if (fairight[i] < 1e-7)
		{
			alpharight[i] = 0.0;
		}
		else
		{
			alpharight[i] = dright[i];
		}
	}

	double alphar = 0.0;

	for (int i = 0; i < 5; i++) alphar += alpharight[i];

	double omegaright[5];
	for (int i = 0; i < 5; i++) omegaright[i] = alpharight[i] / alphar;

	double right = 0.0;
	for (int i = 0; i < 5; i++) right += omegaright[i] * qright[i];

	var = right;
	der2 = 0.0;
	der1 = 0.0;
}

void TENO7_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)
{

	if (block.uniform == true)
	{
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 2))
		{
			TENO7(left.line.right, right.line.left,
				fluids[-3 * block.ny].convar, fluids[-2 * block.ny].convar, fluids[-block.ny].convar,
				fluids[0].convar,
				fluids[block.ny].convar, fluids[2 * block.ny].convar, fluids[3 * block.ny].convar,
				fluids[0].dx,block);
		}

		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 2))
		{
			double wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4];
			YchangetoX(wn3tmp, fluids[-3].convar);
			YchangetoX(wn2tmp, fluids[-2].convar);
			YchangetoX(wn1tmp, fluids[-1].convar);
			YchangetoX(wtmp, fluids[0].convar);
			YchangetoX(wp1tmp, fluids[1].convar);
			YchangetoX(wp2tmp, fluids[2].convar);
			YchangetoX(wp3tmp, fluids[3].convar);

			TENO7(down.line.right, up.line.left,
				wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp,
				fluids[0].dy, block);
		}
	}
	else
	{
		// for non-uniform mesh.
		Point2d voidpoint;
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 2))
		{
			double dx = fluids[0].dx;
			double normal[2];
			double wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4];

			// cell left reconstruction
			Copy_Array(normal, left.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3 * block.ny].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2 * block.ny].convar, normal);
			Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);
			Global_to_Local(wp2tmp, fluids[2 * block.ny].convar, normal);
			Global_to_Local(wp3tmp, fluids[3 * block.ny].convar, normal);
			TENO7(left.line.right, voidpoint, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dx, block);

			// cell right reconstruction
			Copy_Array(normal, right.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3 * block.ny].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2 * block.ny].convar, normal);
			Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);
			Global_to_Local(wp2tmp, fluids[2 * block.ny].convar, normal);
			Global_to_Local(wp3tmp, fluids[3 * block.ny].convar, normal);
			TENO7(voidpoint, right.line.left, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dx, block);
		}

		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 2))
		{
			double dy = fluids[0].dy;
			double normal[2];
			double wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4];

			// down interface reconstruction
			Copy_Array(normal, down.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2].convar, normal);
			Global_to_Local(wn1tmp, fluids[-1].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[1].convar, normal);
			Global_to_Local(wp2tmp, fluids[2].convar, normal);
			Global_to_Local(wp3tmp, fluids[3].convar, normal);
			TENO7(down.line.right, voidpoint, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dy, block);

			// up interface reconstruction
			Copy_Array(normal, up.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2].convar, normal);
			Global_to_Local(wn1tmp, fluids[-1].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[1].convar, normal);
			Global_to_Local(wp2tmp, fluids[2].convar, normal);
			Global_to_Local(wp3tmp, fluids[3].convar, normal);
			TENO7(voidpoint, up.line.left, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dy, block);
		}
	}
}
void TENO7(Point2d& left, Point2d& right, double* wn3, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double* wp3, double h, Block2d block)
{
	double ren3[4], ren2[4], ren1[4], re0[4], rep1[4], rep2[4], rep3[4];
	double var[4], der1[4], der2[4];

	double base_left[4];
	double base_right[4];
	double wn3_primvar[4], wn2_primvar[4], wn1_primvar[4], w_primvar[4], wp1_primvar[4], wp2_primvar[4], wp3_primvar[4];
	//Convar_to_primvar_2D(wn3_primvar, wn3);
	//Convar_to_primvar_2D(wn2_primvar, wn2);
	Convar_to_primvar_2D(wn1_primvar, wn1);
	Convar_to_primvar_2D(w_primvar, w);
	Convar_to_primvar_2D(wp1_primvar, wp1);
	//Convar_to_primvar_2D(wp2_primvar, wp2);
	//Convar_to_primvar_2D(wp3_primvar, wp3);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
		base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren3[i] = wn3[i];
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
			rep3[i] = wp3[i];
		}
	}
	else
	{
		Convar_to_char(ren3, base_left, wn3);
		Convar_to_char(ren2, base_left, wn2);
		Convar_to_char(ren1, base_left, wn1);
		Convar_to_char(re0, base_left, w);
		Convar_to_char(rep1, base_left, wp1);
		Convar_to_char(rep2, base_left, wp2);
		Convar_to_char(rep3, base_left, wp3);
	}


	if (left.gama > 1e-20) // >0
	{
		for (int i = 0; i < 4; i++)
		{
			left.sensor = 1.0;
			TENO7_left(left.gama,var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], h);
		}
	}
	else
	{
		for (int i = 0; i < 4; i++)
		{
			left.sensor = 0.0;
			TENO7opt_left(var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], h);
		}

	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			left.convar[i] = var[i];
			left.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(left.convar, base_left, var);
		Char_to_convar(left.der1x, base_left, der1);
	}



	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren3[i] = wn3[i];
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
			rep3[i] = wp3[i];
		}
	}
	else
	{
		Convar_to_char(ren3, base_right, wn3);
		Convar_to_char(ren2, base_right, wn2);
		Convar_to_char(ren1, base_right, wn1);
		Convar_to_char(re0, base_right, w);
		Convar_to_char(rep1, base_right, wp1);
		Convar_to_char(rep2, base_right, wp2);
		Convar_to_char(rep3, base_right, wp3);
	}

	if (right.gama > 1e-20)
	{
		for (int i = 0; i < 4; i++)
		{
			right.sensor = 1.0;
			TENO7_right(right.gama,var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], h);
		}
	}
	else
	{
		for (int i = 0; i < 4; i++)
		{
			right.sensor = 0.0;
			TENO7opt_right(var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], h);
		}

	}
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			right.convar[i] = var[i];
			right.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(right.convar, base_right, var);
		Char_to_convar(right.der1x, base_right, der1);
	}

	Check_Order_Reduce_by_Lambda_2D(right.is_reduce_order, right.convar);
	Check_Order_Reduce_by_Lambda_2D(left.is_reduce_order, left.convar);

	if (left.is_reduce_order == true || right.is_reduce_order == true)
	{
		if (is_reduce_order_warning == true)
			cout << " TENO6-cell-splitting order reduce" << endl;
		for (int m = 0; m < 4; m++)
		{
			right.convar[m] = w[m];
			left.convar[m] = w[m];
			//right.der1x[m] = 0.0;
			//left.der1x[m] = 0.0;
		}
	}
}
void TENO7_left(double& gama, double& var, double& der1, double& der2, double wn3, double wn2, double wn1, double w, double wp1, double wp2, double wp3, double h)
{
	double qleft[5];
	double dleft[5] = { 18.0 / 35.0, 3.0 / 35.0, 9.0 / 35.0, 1.0 / 35.0 , 4.0 / 35.0 };

	qleft[0] = (2.0 / 6.0) * wn1 + (5.0 / 6.0) * w + (-1.0 / 6.0) * wp1;
	qleft[1] = (11.0 / 6.0) * w + (-7.0 / 6.0) * wp1 + (2.0 / 6.0) * wp2;
	qleft[2] = (-1.0 / 6.0) * wn2 + (5.0 / 6.0) * wn1 + (2.0 / 6.0) * w;
	qleft[3] = (25.0 / 12.0) * w + (-23.0 / 12.0) * wp1 + (13.0 / 12.0) * wp2 + (-3.0 / 12.0) * wp3;
	qleft[4] = (1.0 / 12.0) * wn3 + (-5.0 / 12.0) * wn2 + (13.0 / 12.0) * wn1 + (3.0 / 12.0) * w;


	double beta[5];

	beta[0] = (13.0 / 12.0) * pow(wn1 - 2.0 * w + wp1, 2)
		+ 0.25 * pow(wn1 - wp1, 2);

	beta[1] = (13.0 / 12.0) * pow(w - 2.0 * wp1 + wp2, 2)
		+ 0.25 * pow(3.0 * w - 4.0 * wp1 + wp2, 2);

	beta[2] = (13.0 / 12.0) * pow(wn2 - 2.0 * wn1 + w, 2)
		+ 0.25 * pow(wn2 - 4.0 * wn1 + 3.0 * w, 2);

	beta[3] = (2107.0 * w * w - 9402.0 * w * wp1 + 7042.0 * w * wp2 - 1854.0 * w * wp3
		+ 11003.0 * wp1 * wp1 - 17246.0 * wp1 * wp2 + 4642.0 * wp1 * wp3
		+ 7043.0 * wp2 * wp2 - 3882.0 * wp2 * wp3
		+ 547.0 * wp3 * wp3) / 240.0;

	beta[4] = (547.0 * wn3 * wn3 - 3882.0 * wn3 * wn2 + 4642.0 * wn3 * wn1 - 1854.0 * wn3 * w
		+ 7043.0 * wn2 * wn2 - 17246.0 * wn2 * wn1 + 7042.0 * wn2 * w
		+ 11003.0 * wn1 * wn1 - 9402.0 * wn1 * w
		+ 2107.0 * w * w) / 240.0;


	double epsilon2 = 1e-40;

	double alphaleft[5];

	for (int i = 0; i < 5; i++)
	{
		if (beta[i] < 1e-30) beta[i] = 0.0;
	}


	double gamaleft[5];
	double tau7 = fabs(beta[4] - beta[3]);


	for (int i = 0; i < 5; i++)
	{
		gamaleft[i] = pow((1.0 + tau7 / (beta[i] + epsilon2)), 6);
	}
	double gamal = 0.0;
	for (int i = 0; i < 5; i++)
	{
		gamal += gamaleft[i];
	}
	double faileft[5];
	for (int i = 0; i < 5; i++)
	{
		faileft[i] = gamaleft[i] / gamal;
	}


	double m = 8.0 - floor(gama * 1.0);
	double CT = pow(10.0, -m);
	//double CT = pow(10.0, -7);
	//cout << gama << endl;

	for (int i = 0; i < 5; i++)
	{
		if (faileft[i] < CT)
		{
			alphaleft[i] = 0.0;
		}
		else
		{
			alphaleft[i] = dleft[i];
		}
	}


	double alphal = 0.0;
	for (int i = 0; i < 5; i++) alphal += alphaleft[i];

	double omegaleft[5];
	for (int i = 0; i < 5; i++) omegaleft[i] = alphaleft[i] / alphal;

	double left = 0.0;
	for (int i = 0; i < 5; i++) left += omegaleft[i] * qleft[i];


	var = left;
	der2 = 0.0;
	der1 = 0.0;
}
void TENO7_right(double& gama,double& var, double& der1, double& der2, double wn3, double wn2, double wn1, double w, double wp1, double wp2, double wp3, double h)
{
	double qright[5];
	double dright[5] = { 18.0 / 35.0, 9.0 / 35.0, 3.0 / 35.0, 4.0 / 35.0 , 1.0 / 35.0 };

	qright[0] = (-1.0 / 6.0) * wn1 + (5.0 / 6.0) * w + (2.0 / 6.0) * wp1;
	qright[1] = (2.0 / 6.0) * w + (5.0 / 6.0) * wp1 - (1.0 / 6.0) * wp2;
	qright[2] = (2.0 / 6.0) * wn2 - (7.0 / 6.0) * wn1 + (11.0 / 6.0) * w;
	qright[3] = (3.0 / 12.0) * w + (13.0 / 12.0) * wp1 - (5.0 / 12.0) * wp2 + (1.0 / 12.0) * wp3;
	qright[4] = (-3.0 / 12.0) * wn3 + (13.0 / 12.0) * wn2 - (23.0 / 12.0) * wn1 + (25.0 / 12.0) * w;


	double beta[5];

	beta[0] = (13.0 / 12.0) * pow(wn1 - 2.0 * w + wp1, 2)
		+ 0.25 * pow(wn1 - wp1, 2);

	beta[1] = (13.0 / 12.0) * pow(w - 2.0 * wp1 + wp2, 2)
		+ 0.25 * pow(3.0 * w - 4.0 * wp1 + wp2, 2);

	beta[2] = (13.0 / 12.0) * pow(wn2 - 2.0 * wn1 + w, 2)
		+ 0.25 * pow(wn2 - 4.0 * wn1 + 3.0 * w, 2);


	beta[3] = (1.0 / 36.0) * pow(-11.0 * w + 18.0 * wp1 - 9.0 * wp2 + 2 * wp3, 2) + (13.0 / 12.0) * pow(2.0 * w - 5.0 * wp1 + 4.0 * wp2 - wp3, 2)
		+ (781.0 / 720.0) * pow(-w + 3.0 * wp1 - 3.0 * wp2 + wp3, 2);


	beta[4] = (547.0 * wn3 * wn3 - 3882.0 * wn3 * wn2 + 4642.0 * wn3 * wn1 - 1854.0 * wn3 * w
		+ 7043.0 * wn2 * wn2 - 17246.0 * wn2 * wn1 + 7042.0 * wn2 * w
		+ 11003.0 * wn1 * wn1 - 9402.0 * wn1 * w
		+ 2107.0 * w * w) / 240.0;



	for (int i = 0; i < 5; i++)
	{
		if (beta[i] < 1e-30) beta[i] = 0.0;
	}

	double epsilon2 = 1e-40;


	double alpharight[5];
	double gamaright[5];

	double tau7 = fabs(beta[4] - beta[3]);

	for (int i = 0; i < 5; i++)
	{
		gamaright[i] = pow((1.0 + tau7 / (beta[i] + epsilon2)), 6);
	}
	double gamar = 0.0;
	for (int i = 0; i < 5; i++)
	{
		gamar += gamaright[i];
	}
	double fairight[5];
	for (int i = 0; i < 5; i++)
	{
		fairight[i] = gamaright[i] / gamar;
	}


	double m = 8.0 - floor(gama * 1.0);
	double CT = pow(10.0, -m);
	//double CT = pow(10.0, -7);
	//cout << gama << endl;

	for (int i = 0; i < 5; i++)
	{
		if (fairight[i] < CT)
		{
			alpharight[i] = 0.0;
		}
		else
		{
			alpharight[i] = dright[i];
		}
	}

	double alphar = 0.0;

	for (int i = 0; i < 5; i++) alphar += alpharight[i];

	double omegaright[5];
	for (int i = 0; i < 5; i++) omegaright[i] = alpharight[i] / alphar;

	double right = 0.0;
	for (int i = 0; i < 5; i++) right += omegaright[i] * qright[i];

	var = right;
	der2 = 0.0;
	der1 = 0.0;
}

void TENO70_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)
{

	if (block.uniform == true)
	{
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 2))
		{
			TENO70(left.line.right, right.line.left,
				fluids[-3 * block.ny].convar, fluids[-2 * block.ny].convar, fluids[-block.ny].convar,
				fluids[0].convar,
				fluids[block.ny].convar, fluids[2 * block.ny].convar, fluids[3 * block.ny].convar,
				fluids[0].dx);
		}

		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 2))
		{
			double wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4];
			YchangetoX(wn3tmp, fluids[-3].convar);
			YchangetoX(wn2tmp, fluids[-2].convar);
			YchangetoX(wn1tmp, fluids[-1].convar);
			YchangetoX(wtmp, fluids[0].convar);
			YchangetoX(wp1tmp, fluids[1].convar);
			YchangetoX(wp2tmp, fluids[2].convar);
			YchangetoX(wp3tmp, fluids[3].convar);

			TENO70(down.line.right, up.line.left,
				wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp,
				fluids[0].dy);
		}
	}
	else
	{
		// for non-uniform mesh.
		Point2d voidpoint;
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 2))
		{
			double dx = fluids[0].dx;
			double normal[2];
			double wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4];

			// cell left reconstruction
			Copy_Array(normal, left.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3 * block.ny].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2 * block.ny].convar, normal);
			Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);
			Global_to_Local(wp2tmp, fluids[2 * block.ny].convar, normal);
			Global_to_Local(wp3tmp, fluids[3 * block.ny].convar, normal);
			TENO70(left.line.right, voidpoint, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dx);

			// cell right reconstruction
			Copy_Array(normal, right.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3 * block.ny].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2 * block.ny].convar, normal);
			Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);
			Global_to_Local(wp2tmp, fluids[2 * block.ny].convar, normal);
			Global_to_Local(wp3tmp, fluids[3 * block.ny].convar, normal);
			TENO70(voidpoint, right.line.left, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dx);
		}

		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 2))
		{
			double dy = fluids[0].dy;
			double normal[2];
			double wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4];

			// down interface reconstruction
			Copy_Array(normal, down.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2].convar, normal);
			Global_to_Local(wn1tmp, fluids[-1].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[1].convar, normal);
			Global_to_Local(wp2tmp, fluids[2].convar, normal);
			Global_to_Local(wp3tmp, fluids[3].convar, normal);
			TENO70(down.line.right, voidpoint, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dy);

			// up interface reconstruction
			Copy_Array(normal, up.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2].convar, normal);
			Global_to_Local(wn1tmp, fluids[-1].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[1].convar, normal);
			Global_to_Local(wp2tmp, fluids[2].convar, normal);
			Global_to_Local(wp3tmp, fluids[3].convar, normal);
			TENO70(voidpoint, up.line.left, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dy);
		}
	}


}
void TENO70(Point2d& left, Point2d& right, double* wn3, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double* wp3, double h)
{
	double ren3[4], ren2[4], ren1[4], re0[4], rep1[4], rep2[4], rep3[4];
	double var[4], der1[4], der2[4];

	double base_left[4];
	double base_right[4];
	double wn3_primvar[4], wn2_primvar[4], wn1_primvar[4], w_primvar[4], wp1_primvar[4], wp2_primvar[4], wp3_primvar[4];
	//Convar_to_primvar_2D(wn3_primvar, wn3);
	//Convar_to_primvar_2D(wn2_primvar, wn2);
	Convar_to_primvar_2D(wn1_primvar, wn1);
	Convar_to_primvar_2D(w_primvar, w);
	Convar_to_primvar_2D(wp1_primvar, wp1);
	//Convar_to_primvar_2D(wp2_primvar, wp2);
	//Convar_to_primvar_2D(wp3_primvar, wp3);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
		base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren3[i] = wn3[i];
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
			rep3[i] = wp3[i];
		}
	}
	else
	{
		Convar_to_char(ren3, base_left, wn3);
		Convar_to_char(ren2, base_left, wn2);
		Convar_to_char(ren1, base_left, wn1);
		Convar_to_char(re0, base_left, w);
		Convar_to_char(rep1, base_left, wp1);
		Convar_to_char(rep2, base_left, wp2);
		Convar_to_char(rep3, base_left, wp3);
	}


	for (int i = 0; i < 4; i++)
	{
		TENO70_left(var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			left.convar[i] = var[i];
			left.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(left.convar, base_left, var);
		for (int i = 0; i < 4; i++)
		{
			left.der1x[i] = der1[i];
		}
	}



	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren3[i] = wn3[i];
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
			rep3[i] = wp3[i];
		}
	}
	else
	{
		Convar_to_char(ren3, base_right, wn3);
		Convar_to_char(ren2, base_right, wn2);
		Convar_to_char(ren1, base_right, wn1);
		Convar_to_char(re0, base_right, w);
		Convar_to_char(rep1, base_right, wp1);
		Convar_to_char(rep2, base_right, wp2);
		Convar_to_char(rep3, base_right, wp3);
	}


	for (int i = 0; i < 4; i++)
	{
		TENO70_right(var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			right.convar[i] = var[i];
			right.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(right.convar, base_right, var);
		for (int i = 0; i < 4; i++)
		{
			right.der1x[i] = der1[i];
		}
	}

	Check_Order_Reduce_by_Lambda_2D(right.is_reduce_order, right.convar);
	Check_Order_Reduce_by_Lambda_2D(left.is_reduce_order, left.convar);

	if (left.is_reduce_order == true || right.is_reduce_order == true)
	{
		if (is_reduce_order_warning == true)
			cout << " TENO6-cell-splitting order reduce" << endl;
		for (int m = 0; m < 4; m++)
		{
			right.convar[m] = w[m];
			left.convar[m] = w[m];
		}
	}
}
void TENO70_left(double& var, double& der1, double& der2, double wn3, double wn2, double wn1, double w, double wp1, double wp2, double wp3, double h)
{
	double qleft[5];
	double dleft[5] = { 18.0 / 35.0, 3.0 / 35.0, 9.0 / 35.0, 1.0 / 35.0 , 4.0 / 35.0 };

	qleft[0] = (2.0 / 6.0) * wn1 + (5.0 / 6.0) * w + (-1.0 / 6.0) * wp1;
	qleft[1] = (11.0 / 6.0) * w + (-7.0 / 6.0) * wp1 + (2.0 / 6.0) * wp2;
	qleft[2] = (-1.0 / 6.0) * wn2 + (5.0 / 6.0) * wn1 + (2.0 / 6.0) * w;
	qleft[3] = (25.0 / 12.0) * w + (-23.0 / 12.0) * wp1 + (13.0 / 12.0) * wp2 + (-3.0 / 12.0) * wp3;
	qleft[4] = (1.0 / 12.0) * wn3 + (-5.0 / 12.0) * wn2 + (13.0 / 12.0) * wn1 + (3.0 / 12.0) * w;


	double beta[5];

	beta[0] = (13.0 / 12.0) * pow(wn1 - 2.0 * w + wp1, 2)
		+ 0.25 * pow(wn1 - wp1, 2);

	beta[1] = (13.0 / 12.0) * pow(w - 2.0 * wp1 + wp2, 2)
		+ 0.25 * pow(3.0 * w - 4.0 * wp1 + wp2, 2);

	beta[2] = (13.0 / 12.0) * pow(wn2 - 2.0 * wn1 + w, 2)
		+ 0.25 * pow(wn2 - 4.0 * wn1 + 3.0 * w, 2);

	beta[3] = (2107.0 * w * w - 9402.0 * w * wp1 + 7042.0 * w * wp2 - 1854.0 * w * wp3
		+ 11003.0 * wp1 * wp1 - 17246.0 * wp1 * wp2 + 4642.0 * wp1 * wp3
		+ 7043.0 * wp2 * wp2 - 3882.0 * wp2 * wp3
		+ 547.0 * wp3 * wp3) / 240.0;

	beta[4] = (547.0 * wn3 * wn3 - 3882.0 * wn3 * wn2 + 4642.0 * wn3 * wn1 - 1854.0 * wn3 * w
		+ 7043.0 * wn2 * wn2 - 17246.0 * wn2 * wn1 + 7042.0 * wn2 * w
		+ 11003.0 * wn1 * wn1 - 9402.0 * wn1 * w
		+ 2107.0 * w * w) / 240.0;


	double epsilon2 = 1e-40;

	double alphaleft[5];


	for (int i = 0; i < 5; i++)
	{
		if (beta[i] < 1e-30) beta[i] = 0.0;
	}


	double gamaleft[5];
	double tau7 = fabs(beta[4] - beta[3]);

	for (int i = 0; i < 5; i++)
	{
		gamaleft[i] = pow((1.0 + tau7 / (beta[i] + epsilon2)), 6);
	}
	double gamal = 0.0;
	for (int i = 0; i < 5; i++)
	{
		gamal += gamaleft[i];
	}
	double faileft[5];
	for (int i = 0; i < 5; i++)
	{
		faileft[i] = gamaleft[i] / gamal;
	}
	for (int i = 0; i < 5; i++)
	{
		if (faileft[i] < 1e-7)
		{
			alphaleft[i] = 0.0;
		}
		else
		{
			alphaleft[i] = dleft[i];
		}
	}


	double alphal = 0.0;
	for (int i = 0; i < 5; i++) alphal += alphaleft[i];

	double omegaleft[5];
	for (int i = 0; i < 5; i++) omegaleft[i] = alphaleft[i] / alphal;

	double left = 0.0;
	for (int i = 0; i < 5; i++) left += omegaleft[i] * qleft[i];


	var = left;
	der2 = 0.0;

	double aa = fabs((omegaleft[0] / dleft[0]) - 1.0) + fabs((omegaleft[1] / dleft[1]) - 1.0) + fabs((omegaleft[2] / dleft[2]) - 1.0) + fabs((omegaleft[3] / dleft[3]) - 1.0) + fabs((omegaleft[4] / dleft[4]) - 1.0);
	double mind = min({ dleft[0], dleft[1], dleft[2] , dleft[3] , dleft[4] });
	double bb = fabs((1.0 / mind) - 1.0) + 4.0;
	der1 = aa / bb;

	if (tau7 > 1e-7)
	{
		double aa = fabs((omegaleft[0] / dleft[0]) - 1.0) + fabs((omegaleft[1] / dleft[1]) - 1.0) + fabs((omegaleft[2] / dleft[2]) - 1.0) + fabs((omegaleft[3] / dleft[3]) - 1.0) + fabs((omegaleft[4] / dleft[4]) - 1.0);
		double mind = min({ dleft[0], dleft[1], dleft[2] , dleft[3] , dleft[4] });
		double bb = fabs((1.0 / mind) - 1.0) + 4.0;
		der1 = aa / bb;
	}
	else
	{
		der1 = 0.0;
	}
}
void TENO70_right(double& var, double& der1, double& der2, double wn3, double wn2, double wn1, double w, double wp1, double wp2, double wp3, double h)
{
	double qright[5];
	double dright[5] = { 18.0 / 35.0, 9.0 / 35.0, 3.0 / 35.0, 4.0 / 35.0 , 1.0 / 35.0 };

	qright[0] = (-1.0 / 6.0) * wn1 + (5.0 / 6.0) * w + (2.0 / 6.0) * wp1;
	qright[1] = (2.0 / 6.0) * w + (5.0 / 6.0) * wp1 - (1.0 / 6.0) * wp2;
	qright[2] = (2.0 / 6.0) * wn2 - (7.0 / 6.0) * wn1 + (11.0 / 6.0) * w;
	qright[3] = (3.0 / 12.0) * w + (13.0 / 12.0) * wp1 - (5.0 / 12.0) * wp2 + (1.0 / 12.0) * wp3;
	qright[4] = (-3.0 / 12.0) * wn3 + (13.0 / 12.0) * wn2 - (23.0 / 12.0) * wn1 + (25.0 / 12.0) * w;


	double beta[5];

	beta[0] = (13.0 / 12.0) * pow(wn1 - 2.0 * w + wp1, 2)
		+ 0.25 * pow(wn1 - wp1, 2);

	beta[1] = (13.0 / 12.0) * pow(w - 2.0 * wp1 + wp2, 2)
		+ 0.25 * pow(3.0 * w - 4.0 * wp1 + wp2, 2);

	beta[2] = (13.0 / 12.0) * pow(wn2 - 2.0 * wn1 + w, 2)
		+ 0.25 * pow(wn2 - 4.0 * wn1 + 3.0 * w, 2);


	beta[3] = (1.0 / 36.0) * pow(-11.0 * w + 18.0 * wp1 - 9.0 * wp2 + 2 * wp3, 2) + (13.0 / 12.0) * pow(2.0 * w - 5.0 * wp1 + 4.0 * wp2 - wp3, 2)
		+ (781.0 / 720.0) * pow(-w + 3.0 * wp1 - 3.0 * wp2 + wp3, 2);


	beta[4] = (547.0 * wn3 * wn3 - 3882.0 * wn3 * wn2 + 4642.0 * wn3 * wn1 - 1854.0 * wn3 * w
		+ 7043.0 * wn2 * wn2 - 17246.0 * wn2 * wn1 + 7042.0 * wn2 * w
		+ 11003.0 * wn1 * wn1 - 9402.0 * wn1 * w
		+ 2107.0 * w * w) / 240.0;



	for (int i = 0; i < 5; i++)
	{
		if (beta[i] < 1e-30) beta[i] = 0.0;
	}

	double epsilon2 = 1e-40;


	double alpharight[5];
	double gamaright[5];

	double tau7 = fabs(beta[4] - beta[3]);

	for (int i = 0; i < 5; i++)
	{
		gamaright[i] = pow((1.0 + tau7 / (beta[i] + epsilon2)), 6);
	}
	double gamar = 0.0;
	for (int i = 0; i < 5; i++)
	{
		gamar += gamaright[i];
	}
	double fairight[5];
	for (int i = 0; i < 5; i++)
	{
		fairight[i] = gamaright[i] / gamar;
	}
	for (int i = 0; i < 5; i++)
	{
		if (fairight[i] < 1e-7)
		{
			alpharight[i] = 0.0;
		}
		else
		{
			alpharight[i] = dright[i];
		}
	}

	double alphar = 0.0;

	for (int i = 0; i < 5; i++) alphar += alpharight[i];

	double omegaright[5];
	for (int i = 0; i < 5; i++) omegaright[i] = alpharight[i] / alphar;

	double right = 0.0;
	for (int i = 0; i < 5; i++) right += omegaright[i] * qright[i];

	var = right;
	der2 = 0.0;

	if (tau7 > 1e-7)
	{
		double aa = fabs((omegaright[0] / dright[0]) - 1.0) + fabs((omegaright[1] / dright[1]) - 1.0) + fabs((omegaright[2] / dright[2]) - 1.0) + fabs((omegaright[3] / dright[3]) - 1.0) + fabs((omegaright[4] / dright[4]) - 1.0);
		double mind = min({ dright[0], dright[1], dright[2] , dright[3] , dright[4] });
		double bb = fabs((1.0 / mind) - 1.0) + 4.0;
		der1 = aa / bb;
	}
	else
	{
		der1 = 0.0;
	}
}

void TENO7opt_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)
{

	if (block.uniform == true)
	{
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 2))
		{
			TENO7opt(left.line.right, right.line.left,
				fluids[-3 * block.ny].convar, fluids[-2 * block.ny].convar, fluids[-block.ny].convar,
				fluids[0].convar,
				fluids[block.ny].convar, fluids[2 * block.ny].convar, fluids[3 * block.ny].convar,
				fluids[0].dx);
		}

		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 2))
		{
			double wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4];
			YchangetoX(wn3tmp, fluids[-3].convar);
			YchangetoX(wn2tmp, fluids[-2].convar);
			YchangetoX(wn1tmp, fluids[-1].convar);
			YchangetoX(wtmp, fluids[0].convar);
			YchangetoX(wp1tmp, fluids[1].convar);
			YchangetoX(wp2tmp, fluids[2].convar);
			YchangetoX(wp3tmp, fluids[3].convar);

			TENO7opt(down.line.right, up.line.left,
				wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp,
				fluids[0].dy);
		}
	}
	else
	{
		// for non-uniform mesh.
		Point2d voidpoint;
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 2))
		{
			double dx = fluids[0].dx;
			double normal[2];
			double wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4];

			// cell left reconstruction
			Copy_Array(normal, left.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3 * block.ny].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2 * block.ny].convar, normal);
			Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);
			Global_to_Local(wp2tmp, fluids[2 * block.ny].convar, normal);
			Global_to_Local(wp3tmp, fluids[3 * block.ny].convar, normal);
			TENO7opt(left.line.right, voidpoint, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dx);

			// cell right reconstruction
			Copy_Array(normal, right.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3 * block.ny].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2 * block.ny].convar, normal);
			Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);
			Global_to_Local(wp2tmp, fluids[2 * block.ny].convar, normal);
			Global_to_Local(wp3tmp, fluids[3 * block.ny].convar, normal);
			TENO7opt(voidpoint, right.line.left, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dx);
		}

		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 2))
		{
			double dy = fluids[0].dy;
			double normal[2];
			double wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4];

			// down interface reconstruction
			Copy_Array(normal, down.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2].convar, normal);
			Global_to_Local(wn1tmp, fluids[-1].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[1].convar, normal);
			Global_to_Local(wp2tmp, fluids[2].convar, normal);
			Global_to_Local(wp3tmp, fluids[3].convar, normal);
			TENO7opt(down.line.right, voidpoint, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dy);

			// up interface reconstruction
			Copy_Array(normal, up.normal, 2);
			Global_to_Local(wn3tmp, fluids[-3].convar, normal);
			Global_to_Local(wn2tmp, fluids[-2].convar, normal);
			Global_to_Local(wn1tmp, fluids[-1].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[1].convar, normal);
			Global_to_Local(wp2tmp, fluids[2].convar, normal);
			Global_to_Local(wp3tmp, fluids[3].convar, normal);
			TENO7opt(voidpoint, up.line.left, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, dy);
		}
	}
}
void TENO7opt(Point2d& left, Point2d& right, double* wn3, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double* wp3, double h)
{
	double ren3[4], ren2[4], ren1[4], re0[4], rep1[4], rep2[4], rep3[4];
	double var[4], der1[4], der2[4];

	double base_left[4];
	double base_right[4];
	double wn3_primvar[4], wn2_primvar[4], wn1_primvar[4], w_primvar[4], wp1_primvar[4], wp2_primvar[4], wp3_primvar[4];
	//Convar_to_primvar_2D(wn3_primvar, wn3);
	//Convar_to_primvar_2D(wn2_primvar, wn2);
	Convar_to_primvar_2D(wn1_primvar, wn1);
	Convar_to_primvar_2D(w_primvar, w);
	Convar_to_primvar_2D(wp1_primvar, wp1);
	//Convar_to_primvar_2D(wp2_primvar, wp2);
	//Convar_to_primvar_2D(wp3_primvar, wp3);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
		base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren3[i] = wn3[i];
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
			rep3[i] = wp3[i];
		}
	}
	else
	{
		Convar_to_char(ren3, base_left, wn3);
		Convar_to_char(ren2, base_left, wn2);
		Convar_to_char(ren1, base_left, wn1);
		Convar_to_char(re0, base_left, w);
		Convar_to_char(rep1, base_left, wp1);
		Convar_to_char(rep2, base_left, wp2);
		Convar_to_char(rep3, base_left, wp3);
	}


	for (int i = 0; i < 4; i++)
	{
		TENO7opt_left(var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			left.convar[i] = var[i];
			left.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(left.convar, base_left, var);
		Char_to_convar(left.der1x, base_left, der1);
	}



	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren3[i] = wn3[i];
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
			rep3[i] = wp3[i];
		}
	}
	else
	{
		Convar_to_char(ren3, base_right, wn3);
		Convar_to_char(ren2, base_right, wn2);
		Convar_to_char(ren1, base_right, wn1);
		Convar_to_char(re0, base_right, w);
		Convar_to_char(rep1, base_right, wp1);
		Convar_to_char(rep2, base_right, wp2);
		Convar_to_char(rep3, base_right, wp3);
	}


	for (int i = 0; i < 4; i++)
	{
		TENO7opt_right(var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			right.convar[i] = var[i];
			right.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(right.convar, base_right, var);
		Char_to_convar(right.der1x, base_right, der1);
	}

	Check_Order_Reduce_by_Lambda_2D(right.is_reduce_order, right.convar);
	Check_Order_Reduce_by_Lambda_2D(left.is_reduce_order, left.convar);

	if (left.is_reduce_order == true || right.is_reduce_order == true)
	{
		if (is_reduce_order_warning == true)
			cout << " TENO6-cell-splitting order reduce" << endl;
		for (int m = 0; m < 4; m++)
		{
			right.convar[m] = w[m];
			left.convar[m] = w[m];
			right.der1x[m] = 0.0;
			left.der1x[m] = 0.0;
		}
	}
}
void TENO7opt_left(double& var, double& der1, double& der2, double wn3, double wn2, double wn1, double w, double wp1, double wp2, double wp3, double h)
{
	double qleft[5];

	qleft[0] = (2.0 / 6.0) * wn1 + (5.0 / 6.0) * w + (-1.0 / 6.0) * wp1;
	qleft[1] = (11.0 / 6.0) * w + (-7.0 / 6.0) * wp1 + (2.0 / 6.0) * wp2;
	qleft[2] = (-1.0 / 6.0) * wn2 + (5.0 / 6.0) * wn1 + (2.0 / 6.0) * w;
	qleft[3] = (25.0 / 12.0) * w + (-23.0 / 12.0) * wp1 + (13.0 / 12.0) * wp2 + (-3.0 / 12.0) * wp3;
	qleft[4] = (1.0 / 12.0) * wn3 + (-5.0 / 12.0) * wn2 + (13.0 / 12.0) * wn1 + (3.0 / 12.0) * w;

	var = (18.0 / 35.0) * qleft[0] + (3.0 / 35.0) * qleft[1] + (9.0 / 35.0) * qleft[2] + (1.0 / 35.0) * qleft[3] + (4.0 / 35.0) * qleft[4];

	der2 = 0.0;
	der1 = 0.0;
}
void TENO7opt_right(double& var, double& der1, double& der2, double wn3, double wn2, double wn1, double w, double wp1, double wp2, double wp3, double h)
{
	double qright[5];

	qright[0] = (-1.0 / 6.0) * wn1 + (5.0 / 6.0) * w + (2.0 / 6.0) * wp1;
	qright[1] = (2.0 / 6.0) * w + (5.0 / 6.0) * wp1 - (1.0 / 6.0) * wp2;
	qright[2] = (2.0 / 6.0) * wn2 - (7.0 / 6.0) * wn1 + (11.0 / 6.0) * w;
	qright[3] = (3.0 / 12.0) * w + (13.0 / 12.0) * wp1 - (5.0 / 12.0) * wp2 + (1.0 / 12.0) * wp3;
	qright[4] = (-3.0 / 12.0) * wn3 + (13.0 / 12.0) * wn2 - (23.0 / 12.0) * wn1 + (25.0 / 12.0) * w;

	var = (18.0 / 35.0) * qright[0] + (9.0 / 35.0) * qright[1] + (3.0 / 35.0) * qright[2] + (4.0 / 35.0) * qright[3] + (1.0 / 35.0) * qright[4];

	der2 = 0.0;
	der1 = 0.0;
}

void Update_dec(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage, Interface2d* xinterfaces, Interface2d* yinterfaces)
{

	if (stage > block.stages)
	{
		cout << "wrong middle stage,pls check the time marching setting" << endl;
		exit(0);
	}

	double dt = block.dt;

	Update_with_gauss_dec(fluids, xfluxes, yfluxes, block, stage);
}
void Update_with_gauss_dec(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage)
{
	//Note : calculate the final flux of the cell, in fluids, by the obtained final flux of interface, in xfluxes and yfluxes
#pragma omp parallel for collapse(2)
	for (int i = block.ghost - 2; i < block.nodex + block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.nodey + block.ghost + 2; j++)
		{
			int face = i * (block.ny + 1) + j;
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				Local_to_Global(yfluxes[face][stage].gauss[num_gauss].f, yfluxes[face][stage].gauss[num_gauss].normal);
				Local_to_Global(xfluxes[face][stage].gauss[num_gauss].f, xfluxes[face][stage].gauss[num_gauss].normal);
			}
		}
	}

#pragma omp parallel for collapse(2)
	for (int i = block.ghost - 2; i < block.nodex + block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.nodey + block.ghost + 2; j++)
		{
			int cell = i * (block.ny) + j;
			int face = i * (block.ny + 1) + j;

			for (int var = 0; var < 4; var++)
			{
				double total_flux = 0.0;
				for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
				{
					total_flux += gauss_weight[num_gauss] * yfluxes[face][stage].gauss[num_gauss].length * yfluxes[face][stage].gauss[num_gauss].f[var];
					total_flux += -gauss_weight[num_gauss] * yfluxes[face + 1][stage].gauss[num_gauss].length * yfluxes[face + 1][stage].gauss[num_gauss].f[var];
					total_flux += gauss_weight[num_gauss] * xfluxes[face][stage].gauss[num_gauss].length * xfluxes[face][stage].gauss[num_gauss].f[var];
					total_flux += -gauss_weight[num_gauss] * xfluxes[face + block.ny + 1][stage].gauss[num_gauss].length * xfluxes[face + block.ny + 1][stage].gauss[num_gauss].f[var];
				}

				fluids[cell].D_du_ut[var] = -total_flux / fluids[cell].area;
			}
		}
	}
}
void comput_du_dt(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage, Interface2d* xinterfaces, Interface2d* yinterfaces)
{
#pragma omp parallel for collapse(2)
	for (int i = block.ghost - 2; i < block.nodex + block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.nodey + block.ghost + 2; j++)
		{
			int face = i * (block.ny + 1) + j;
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				Local_to_Global(yfluxes[face][stage].gauss[num_gauss].f, yfluxes[face][stage].gauss[num_gauss].normal);
				Local_to_Global(xfluxes[face][stage].gauss[num_gauss].f, xfluxes[face][stage].gauss[num_gauss].normal);
			}
		}
	}

	double dt = block.dt;
	double du, ddu, dux, duy;

#pragma omp parallel for collapse(2)
	for (int i = block.ghost - 2; i < block.nodex + block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.nodey + block.ghost + 2; j++)
		{
			int cell = i * (block.ny) + j;

			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				for (int var = 0; var < 4; var++)
				{

					fluids[cell].D_du_ut[var] = (xfluxes[(i + 1) * (block.ny + 1) + j][stage].gauss[num_gauss].f[var] - xfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].f[var]) / block.dx; // y1 du/dx
				}
			}
		}
	}
}
void Update_RK44(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage, Interface2d* xinterfaces, Interface2d* yinterfaces) // RK33
{
#pragma omp parallel for collapse(2)
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost + 1; j++)
		{
			int face = i * (block.ny + 1) + j;
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				Local_to_Global(yfluxes[face][stage].gauss[num_gauss].f, yfluxes[face][stage].gauss[num_gauss].normal);
				Local_to_Global(xfluxes[face][stage].gauss[num_gauss].f, xfluxes[face][stage].gauss[num_gauss].normal);
			}
		}
	}

	double local;
	double du, un0, un1, du0, un2, dux, duy;
	double dt = block.dt;

#pragma omp parallel for collapse(2)
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			int cell = i * block.ny + j;

			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				for (int var = 0; var < 4; var++)
				{
					double dux = (xfluxes[(i + 1) * (block.ny + 1) + j][stage].gauss[num_gauss].f[var] -
						xfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].f[var]) / block.dx;
					double duy = (yfluxes[i * (block.ny + 1) + j + 1][stage].gauss[num_gauss].f[var] -
						yfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].f[var]) / block.dy;
					double du = dux + duy;

					double local;
					if (stage == 0)
					{
						double un0 = fluids[cell].D_u[var];
						local = un0 - 0.5 * dt * du;  // y1 u

						fluids[cell].D_u[var] = local;
						fluids[cell].D_u0[var] = un0;
						fluids[cell].D_du1[var] = du;

						fluids[cell].convar[var] = local;
						fluids[cell].convar_old[var] = un0;
					}
					else if (stage == 1)
					{
						double un0 = fluids[cell].D_u0[var];
						local = un0 - 0.5 * dt * du;  // y1 u

						fluids[cell].D_u[var] = local;
						fluids[cell].D_u1[var] = un1;
						fluids[cell].D_du2[var] = du;

						fluids[cell].convar[var] = local;
						fluids[cell].convar_old[var] = un0;
					}
					else if (stage == 2)
					{
						double un0 = fluids[cell].D_u0[var];
						local = un0 - dt * du;  // y1 u

						fluids[cell].D_u[var] = local;
						fluids[cell].D_u1[var] = un1;
						fluids[cell].D_du3[var] = du;

						fluids[cell].convar[var] = local;
						fluids[cell].convar_old[var] = un0;
					}
					else // stage == 3
					{
						double un0 = fluids[cell].D_u0[var];
						double du1 = fluids[cell].D_du1[var];
						double du2 = fluids[cell].D_du2[var];
						double du3 = fluids[cell].D_du3[var];
						double du4 = du;
						local = un0 - dt * (1.0 / 6.0 * du1 + 1.0 / 3.0 * du2 + 1.0 / 3.0 * du3 + 1.0 / 6.0 * du4); // n+1 u

						fluids[cell].D_u[var] = local;
						fluids[cell].convar[var] = local;
						fluids[cell].convar_old[var] = un0;
					}
				}
			}
		}
	}
}
void Update_DeC7_1(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, Interface2d* xinterfaces, Interface2d* yinterfaces) // RK22
{
	double local1, local2, local3, local4;
	double un0, F0, F1, F2, F3, F4;
	double dt = block.dt;

	//coefficient
	
	// m = 1
	double theta10 = 0.067728432186156901;
	double theta11 = 0.119744769343411689;
	double theta12 = -0.021735721866558113;
	double theta13 = 0.010635824225415492;
	double theta14 = -0.003700139242414531;

	// m = 2
	double theta20 = 0.040625000000000001;
	double theta21 = 0.303184183323042755;
	double theta22 = 0.177777777777777785;
	double theta23 = -0.030961961100820556;
	double theta24 = 0.009375000000000000;

	// m = 3
	double theta30 = 0.053700139242414534;
	double theta31 = 0.261586397996806719;
	double theta32 = 0.377291277422113658;
	double theta33 = 0.152477452878810538;
	double theta34 = -0.017728432186156898;

	// m = 4
	double theta40 = 0.050000000000000003;
	double theta41 = 0.272222222222222199;
	double theta42 = 0.355555555555555569;
	double theta43 = 0.272222222222222199;
	double theta44 = 0.050000000000000003;


	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			int cell = i * (block.ny) + j;
			for (int var = 0; var < 4; var++)
			{
				un0 = fluids[cell].convar[var];
				F0 = fluids[cell].D_du_ut[var];

				local1 = un0 - dt * F0 * (theta10 + theta11 + theta12 + theta13 + theta14);
				local2 = un0 - dt * F0 * (theta20 + theta21 + theta22 + theta23 + theta24);
				local3 = un0 - dt * F0 * (theta30 + theta31 + theta32 + theta33 + theta34);
				local4 = un0 - dt * F0 * (theta40 + theta41 + theta42 + theta43 + theta44);

				fluids[cell].D_du_ut0[var] = fluids[cell].D_du_ut[var];
				fluids[cell].convar0[var] = fluids[cell].convar[var];
				fluids[cell].convar1[var] = local1;
				fluids[cell].convar2[var] = local2;
				fluids[cell].convar3[var] = local3;
				fluids[cell].convar4[var] = local4; // tn+1
			}
		}
	}
}
void Update_DeC7(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, Interface2d* xinterfaces, Interface2d* yinterfaces) // RK22
{
	double local1, local2, local3, local4;
	double un0, F0, F1, F2, F3, F4;
	double dt = block.dt;

	//coefficient
	
	// m = 1
	double theta10 = 0.067728432186156901;
	double theta11 = 0.119744769343411689;
	double theta12 = -0.021735721866558113;
	double theta13 = 0.010635824225415492;
	double theta14 = -0.003700139242414531;

	// m = 2
	double theta20 = 0.040625000000000001;
	double theta21 = 0.303184183323042755;
	double theta22 = 0.177777777777777785;
	double theta23 = -0.030961961100820556;
	double theta24 = 0.009375000000000000;

	// m = 3
	double theta30 = 0.053700139242414534;
	double theta31 = 0.261586397996806719;
	double theta32 = 0.377291277422113658;
	double theta33 = 0.152477452878810538;
	double theta34 = -0.017728432186156898;

	// m = 4
	double theta40 = 0.050000000000000003;
	double theta41 = 0.272222222222222199;
	double theta42 = 0.355555555555555569;
	double theta43 = 0.272222222222222199;
	double theta44 = 0.050000000000000003;



	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			int cell = i * (block.ny) + j;
			for (int var = 0; var < 4; var++)
			{
				un0 = fluids[cell].convar0[var];
				F0 = fluids[cell].D_du_ut0[var];
				F1 = fluids[cell].D_du_ut1[var];
				F2 = fluids[cell].D_du_ut2[var];
				F3 = fluids[cell].D_du_ut3[var];
				F4 = fluids[cell].D_du_ut4[var];

				local1 = un0 - dt * (theta10 * F0 + theta11 * F1 + theta12 * F2 + theta13 * F3 + theta14 * F4);
				local2 = un0 - dt * (theta20 * F0 + theta21 * F1 + theta22 * F2 + theta23 * F3 + theta24 * F4);
				local3 = un0 - dt * (theta30 * F0 + theta31 * F1 + theta32 * F2 + theta33 * F3 + theta34 * F4);
				local4 = un0 - dt * (theta40 * F0 + theta41 * F1 + theta42 * F2 + theta43 * F3 + theta44 * F4);

				fluids[cell].convar1[var] = local1;
				fluids[cell].convar2[var] = local2;
				fluids[cell].convar3[var] = local3;
				fluids[cell].convar4[var] = local4; // tn+1

			}
		}
	}
}
void Calculate_flux(Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Interface2d* xinterfaces, Interface2d* yinterfaces, Block2d block, int stage, Fluid2d* fluids)
{
#pragma omp parallel for collapse(2)
	for (int i = block.ghost - 2; i < block.nodex + block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.nodey + block.ghost + 2; j++)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				double n1x = fluids[i * block.ny + j].node[0];
				double n1y = fluids[i * block.ny + j].node[1];
				double n2x = fluids[i * block.ny + j + 1].node[0];
				double n2y = fluids[i * block.ny + j + 1].node[1];
				int cell = i * (block.ny + 1) + j;
				flux_function_2d(xfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss], xinterfaces[i * (block.ny + 1) + j].gauss[num_gauss], block.dt, n1x, n1y, n2x, n2y, cell, num_gauss, 0, fluids);
			}
		}
	}

#pragma omp parallel for collapse(2)
	for (int i = block.ghost - 2; i < block.nodex + block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.nodey + block.ghost + 2; j++)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				double n2x = fluids[i * block.ny + j].node[0];
				double n2y = fluids[i * block.ny + j].node[1];
				double n1x = fluids[(i + 1) * block.ny + j].node[0];
				double n1y = fluids[(i + 1) * block.ny + j].node[1];
				int cell = i * (block.ny + 1) + j;
				flux_function_2d(yfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss], yinterfaces[i * (block.ny + 1) + j].gauss[num_gauss], block.dt, n1x, n1y, n2x, n2y, cell, num_gauss, 1, fluids);
			}
		}
	}
}


void Reconstruction_forg0(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{

#pragma omp parallel  for
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			g0reconstruction_2D_normal(&xinterfaces[i * (block.ny + 1) + j], &yinterfaces[i * (block.ny + 1) + j], &fluids[i * (block.ny) + j], block);
		}
	}

	// then get the guass point value. That is, so called multi-dimensional property
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nx - block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost + 1; j++)
		{
			g0reconstruction_2D_tangent(&xinterfaces[i * (block.ny + 1) + j], &yinterfaces[i * (block.ny + 1) + j], &fluids[i * (block.ny) + j], block);
		}
	}
}

void LF2D(Flux2d& flux, Recon2d& interface, double dt)
{
	double pl[4], pr[4];
	Convar_to_primvar_2D(pl, interface.left.convar);
	Convar_to_primvar_2D(pr, interface.right.convar);

	double k[2];
	k[0] = abs(pl[1]) + sqrt(Gamma * pl[3] / pl[0]);
	k[1] = abs(pr[1]) + sqrt(Gamma * pr[3] / pr[0]);
	double beta = k[0];
	if (k[1] > k[0]) { beta = k[1]; }
	double flux_l[4], flux_r[4];
	get_flux(pl, flux_l);
	get_flux(pr, flux_r);

	for (int m = 0; m < 4; m++)
	{
		flux.f[m] = 0.5 * ((flux_l[m] + flux_r[m]) - beta * (interface.right.convar[m] - interface.left.convar[m]));
		flux.f[m] *= dt;
	}
	if (tau_type == NS)
	{
		NS_by_central_difference_prim_2D(flux, interface, dt);
	}

}
void HLLC2D(Flux2d& flux, Recon2d& interface, double dt)
{
	double pl[4], pr[4];
	Convar_to_primvar_2D(pl, interface.left.convar);
	Convar_to_primvar_2D(pr, interface.right.convar);

	double al, ar, pvars, pstar, tmp1, tmp2, tmp3, qk, sl, sr, star;
	al = sqrt(Gamma * pl[3] / pl[0]); //sound speed
	ar = sqrt(Gamma * pr[3] / pr[0]);
	tmp1 = 0.5 * (al + ar);         //avg of sound and density
	tmp2 = 0.5 * (pl[0] + pr[0]);

	pvars = 0.5 * (pl[3] + pr[3]) - 0.5 * (pr[1] - pl[1]) * tmp1 * tmp2;
	pstar = fmax(0.0, pvars);

	double flxtmp[4], qstar[4];

	ESTIME(sl, star, sr, pl[0], pl[1], pl[3], al, pr[0], pr[1], pr[3], ar);

	tmp1 = pr[3] - pl[3] + pl[0] * pl[1] * (sl - pl[1]) - pr[0] * pr[1] * (sr - pr[1]);
	tmp2 = pl[0] * (sl - pl[1]) - pr[0] * (sr - pr[1]);
	star = tmp1 / tmp2;

	if (sl >= 0.0)
	{
		get_flux(pl, flux.f);
	}
	else if (sr <= 0.0)
	{
		get_flux(pr, flux.f);
	}
	else if ((star >= 0.0) && (sl <= 0.0))
	{
		get_flux(pl, flxtmp);
		ustarforHLLC(pl[0], pl[1], pl[2], pl[3], sl, star, qstar);

		for (int m = 0; m < 4; m++)
		{
			flux.f[m] = flxtmp[m] + sl * (qstar[m] - interface.left.convar[m]);
		}
	}
	else if ((star <= 0.0) && (sr >= 0.0))
	{
		get_flux(pr, flxtmp);
		ustarforHLLC(pr[0], pr[1], pr[2], pr[3], sr, star, qstar);
		for (int m = 0; m < 4; m++)
		{
			flux.f[m] = flxtmp[m] + sr * (qstar[m] - interface.right.convar[m]);
		}
	}
	else
	{
		cout << "couldnt be possible that hllc cannot give any result" << endl;
	}

	for (int m = 0; m < 4; m++)
	{
		flux.f[m] *= dt;
	}
	if (tau_type == NS)
	{
		NS_by_central_difference_convar_2D(flux, interface, dt);
	}
}
void NS_by_central_difference_prim_2D(Flux2d& flux, Recon2d& interface, double dt)
{
	// Note : Sixth-order centrial differential scheme for viscous term
	// here convar[i] represents density, u, v, and temperature
	double mu = Mu;
	if (Nu > 0) { mu = Nu * interface.center.convar[0]; }


	Local_to_Global(interface.center.der1x, interface.normal);
	Local_to_Global(interface.center.der1y, interface.normal);
	Local_to_Global(interface.center.convar, interface.normal);

	double u = interface.center.convar[1];
	double v = interface.center.convar[2];

	double ux, uy, vx, vy;
	ux = interface.center.der1x[1];
	vx = interface.center.der1x[2];
	uy = interface.center.der1y[1];
	vy = interface.center.der1y[2];
	double tau_xx = 2 * mu * ux - 2.0 / 3.0 * mu * (ux + vy);
	double tau_xy = mu * (uy + vx);
	double q = u * tau_xx + v * tau_xy + (K + 4) / (2 * Pr) * mu * interface.center.der1x[3];
	flux.f[1] += tau_xx * dt;
	flux.f[2] += tau_xy * dt;
	flux.f[3] += q * dt;
}
void NS_by_central_difference_convar_2D(Flux2d& flux, Recon2d& interface, double dt)
{
	// Note: Write a function calculating the flux of viscous term in NS, used for HLLC Riemann solver
	// If g0type = collisionnless, means the conservative variables used for viscous flux are those of Gauss points, which were interpolated by the line-averaged variables in the tangential direction
	// If g0type = collisionn, means, the conservative variables used for viscous flux are averaged by the left and right variables (line-ageraged)
	// This is because of the first method might cause negative variables for density or pressure, which is not allowed.
	double convar[4];

	//Local_to_Global(interface.left.convar, interface.normal);
	//Local_to_Global(interface.right.convar, interface.normal);
	Local_to_Global(interface.center.der1x, interface.normal);
	Local_to_Global(interface.center.der1y, interface.normal);
	for (int m = 0; m < 4; m++)
	{
		convar[m] = 0.5 * (interface.left.convar[m] + interface.right.convar[m]);
	}

	// here convar[i] represents the conservative variables
	double den = convar[0];
	double mu = Mu;
	if (Nu > 0) { mu = Nu * den; }
	double u = convar[1] / den;
	double v = convar[2] / den;
	double ux, uy, vx, vy;
	ux = (interface.center.der1x[1] - interface.center.der1x[0] * u) / den;
	vx = (interface.center.der1x[2] - interface.center.der1x[0] * v) / den;
	uy = (interface.center.der1y[1] - interface.center.der1y[0] * u) / den;
	vy = (interface.center.der1y[2] - interface.center.der1y[0] * v) / den;
	double tau_xx = 2 * mu * ux - 2.0 / 3.0 * mu * (ux + vy);
	double tau_xy = mu * (uy + vx);
	double Tx = interface.center.der1x[3] / den - u * ux - u * vx
		- convar[3] * interface.center.der1x[0] / den / den;
	double q = u * tau_xx + v * tau_xy + (K + 4) / (2.0 * Pr) * mu * Tx;
	flux.f[1] -= tau_xx * dt;
	flux.f[2] -= tau_xy * dt;
	flux.f[3] -= q * dt;
}

//forward_euler
void S1O1_2D(Block2d& block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.0;

}
void S1O2_2D(Block2d& block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.5;

}
void RK2_2D(Block2d& block)
{
	block.stages = 2;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[1][0][0] = 0.5;
	block.timecoefficient[1][1][0] = 0.5;

}
void S2O4_2D(Block2d& block)
{
	block.stages = 2;
	block.timecoefficient[0][0][0] = 0.5;
	block.timecoefficient[0][0][1] = 1.0 / 8.0;
	block.timecoefficient[1][0][0] = 1.0;
	block.timecoefficient[1][1][0] = 0.0;
	block.timecoefficient[1][0][1] = 1.0 / 6.0;
	block.timecoefficient[1][1][1] = 1.0 / 3.0;

}
void RK4_2D(Block2d& block)
{
	block.stages = 4;
	block.timecoefficient[0][0][0] = 0.5;
	block.timecoefficient[1][1][0] = 0.5;
	block.timecoefficient[2][2][0] = 1.0;
	block.timecoefficient[3][0][0] = 1.0 / 6.0;
	block.timecoefficient[3][1][0] = 1.0 / 3.0;
	block.timecoefficient[3][2][0] = 1.0 / 3.0;
	block.timecoefficient[3][3][0] = 1.0 / 6.0;
}
void Initial_stages(Block2d& block)
{
	for (int i = 0; i < 5; i++) //refers the n stage
	{
		for (int j = 0; j < 5; j++) //refers the nth coefficient at n stage
		{
			for (int k = 0; k < 3; k++) //refers f, derf, der2f
			{
				block.timecoefficient[i][j][k] = 0.0;
			}
		}
	}
	timecoe_list_2d(block);
}



Fluid2d* Setfluid(Block2d& block)
{
	Fluid2d* var = new Fluid2d[block.nx * block.ny]; // dynamic variable (since block.nx is not determined)
	if (var == 0)
	{
		cout << "fluid variable allocate fail...";
		return NULL;
	}
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			var[i * block.ny + j].xindex = i;
			var[i * block.ny + j].yindex = j;
		}
	}

	cout << "fluid variable allocate done..." << endl;
	return var;
}
Interface2d* Setinterface_array(Block2d block)
{
	Interface2d* var = new Interface2d[(block.nx + 1) * (block.ny + 1)];  // dynamic variable (since block.nx is not determined)
	if (var == 0)
	{
		cout << "fluid variable allocate fail...";
		return NULL;
	}
	for (int i = 0; i < block.nx + 1; i++)
	{
		for (int j = 0; j < block.ny + 1; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				var[i * (block.ny + 1) + j].line.left.der1x[k] = 0.0;
				var[i * (block.ny + 1) + j].line.left.der1y[k] = 0.0;

				var[i * (block.ny + 1) + j].line.right.der1x[k] = 0.0;
				var[i * (block.ny + 1) + j].line.right.der1y[k] = 0.0;

				var[i * (block.ny + 1) + j].line.center.der1x[k] = 0.0;
				var[i * (block.ny + 1) + j].line.center.der1y[k] = 0.0;
			}
		}
	}
	cout << "interface variable allocate done..." << endl;
	return var;
}
Flux2d_gauss** Setflux_gauss_array(Block2d block)
{
	Flux2d_gauss** var = new Flux2d_gauss * [(block.nx + 1) * (block.ny + 1)];  // dynamic variable (since block.nx is not determined)

	for (int i = 0; i < block.nx + 1; i++)
	{
		for (int j = 0; j < block.ny + 1; j++)
		{
			// for m th step time marching schemes, m subflux needed
			var[i * (block.ny + 1) + j] = new Flux2d_gauss[block.stages];
		}
	}

	for (int i = 0; i < block.nx + 1; i++)
	{
		for (int j = 0; j < block.ny + 1; j++)
		{
			for (int k = 0; k < block.stages; k++)
			{
				if (gausspoint == 0)
				{
					var[i * (block.ny + 1) + j][k].gauss = new Flux2d[1];
					for (int m = 0; m < 4; m++)
					{
						var[i * (block.ny + 1) + j][k].gauss[0].f[m] = 0.0;
						var[i * (block.ny + 1) + j][k].gauss[0].derf[m] = 0.0;

					}
				}
				else
				{
					var[i * (block.ny + 1) + j][k].gauss = new Flux2d[gausspoint];
					for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
					{
						for (int m = 0; m < 4; m++)
						{
							var[i * (block.ny + 1) + j][k].gauss[num_gauss].f[m] = 0.0;
							var[i * (block.ny + 1) + j][k].gauss[num_gauss].derf[m] = 0.0;
						}
					}
				}
			}
		}
	}
	if (var == 0)
	{
		cout << "fluid variable allocate fail...";
		return NULL;
	}
	cout << "flux with gausspoint variable allocate done..." << endl;
	return var;

}




void SetUniformMesh(Block2d& block, Fluid2d* fluids, Interface2d* xinterfaces, Interface2d* yinterfaces, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes)
{
	//nodex, nodey are the real node
	//interface number = cell number + 1
	block.dx = (block.right - block.left) / block.nodex;
	block.dy = (block.up - block.down) / block.nodey;
	block.overdx = 1 / block.dx;
	block.overdy = 1 / block.dy;

	block.xcell_begin = block.ghost;
	block.xcell_end = block.ghost + block.nodex - 1;
	block.ycell_begin = block.ghost;
	block.ycell_end = block.ghost + block.nodey - 1;

	block.xinterface_begin_n = block.ghost;
	block.xinterface_end_n = block.ghost + block.nodex;
	block.xinterface_begin_t = block.ghost;
	block.xinterface_end_t = block.ghost + block.nodex - 1;

	block.yinterface_begin_n = block.ghost;
	block.yinterface_end_n = block.ghost + block.nodey;
	block.yinterface_begin_t = block.ghost;
	block.yinterface_end_t = block.ghost + block.nodey - 1;

	//cell avg information
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			//two dimension geometry to one dimension store matrix, y direciton first and x direction second
			fluids[i * block.ny + j].dx = block.dx; //cell size
			fluids[i * block.ny + j].dy = block.dy; //cell size
			fluids[i * block.ny + j].coordx = block.left + (i + 0.5 - block.ghost) * block.dx; //cell center location
			fluids[i * block.ny + j].coordy = block.down + (j + 0.5 - block.ghost) * block.dy; //cell center location
			fluids[i * block.ny + j].area = block.dx * block.dy;
			fluids[i * block.ny + j].node[0] = block.left + (i - block.ghost) * block.dx;
			fluids[i * block.ny + j].node[1] = block.down + (j - block.ghost) * block.dy;
			fluids[i * block.ny + j].node[2] = block.left + (i + 1 - block.ghost) * block.dx;
			fluids[i * block.ny + j].node[3] = block.down + (j - block.ghost) * block.dy;
			fluids[i * block.ny + j].node[4] = block.left + (i + 1 - block.ghost) * block.dx;
			fluids[i * block.ny + j].node[5] = block.down + (j + 1 - block.ghost) * block.dy;
			fluids[i * block.ny + j].node[6] = block.left + (i - block.ghost) * block.dx;
			fluids[i * block.ny + j].node[7] = block.down + (j + 1 - block.ghost) * block.dy;
		}
	}

	// interface information
	for (int i = 0; i <= block.nx; i++)
	{
		for (int j = 0; j <= block.ny; j++)
		{
			xinterfaces[i * (block.ny + 1) + j].x = block.left + (i - block.ghost) * block.dx;
			xinterfaces[i * (block.ny + 1) + j].y = block.down + (j - block.ghost + 0.5) * block.dy;
			xinterfaces[i * (block.ny + 1) + j].length = block.dy;
			xinterfaces[i * (block.ny + 1) + j].normal[0] = 1.0;
			xinterfaces[i * (block.ny + 1) + j].normal[1] = 0.0;

			Copy_geo_from_interface_to_line(xinterfaces[i * (block.ny + 1) + j]);
			xinterfaces[i * (block.ny + 1) + j].gauss = new Recon2d[gausspoint];

			Copy_geo_from_interface_to_flux
			(xinterfaces[i * (block.ny + 1) + j], xfluxes[i * (block.ny + 1) + j], block.stages);

			yinterfaces[i * (block.ny + 1) + j].y = block.down + (j - block.ghost) * block.dy;
			yinterfaces[i * (block.ny + 1) + j].x = block.left + (i - block.ghost + 0.5) * block.dx;
			yinterfaces[i * (block.ny + 1) + j].length = block.dx;
			yinterfaces[i * (block.ny + 1) + j].normal[0] = 0.0;
			yinterfaces[i * (block.ny + 1) + j].normal[1] = 1.0;

			Copy_geo_from_interface_to_line(yinterfaces[i * (block.ny + 1) + j]);

			yinterfaces[i * (block.ny + 1) + j].gauss = new Recon2d[gausspoint];

			Copy_geo_from_interface_to_flux
			(yinterfaces[i * (block.ny + 1) + j], yfluxes[i * (block.ny + 1) + j], block.stages);

			Set_Gauss_Coordinate(xinterfaces[i * (block.ny + 1) + j], yinterfaces[i * (block.ny + 1) + j]);
		}
	}
	cout << "set uniform information done..." << endl;
}
void Copy_geo_from_interface_to_line(Interface2d& interface)
{
	interface.line.x = interface.x;
	interface.line.y = interface.y;
	interface.line.normal[0] = interface.normal[0];
	interface.line.normal[1] = interface.normal[1];
}
void Copy_geo_from_interface_to_flux(Interface2d& interface, Flux2d_gauss* flux, int stages)
{
	for (int istage = 0; istage < stages; istage++)
	{
		if (gausspoint == 0)
		{
			int igauss = 0;
			flux[istage].gauss[igauss].normal[0] = interface.normal[0];
			flux[istage].gauss[igauss].normal[1] = interface.normal[1];
			flux[istage].gauss[igauss].length = interface.length;
		}
		else
		{
			for (int igauss = 0; igauss < gausspoint; igauss++)
			{
				flux[istage].gauss[igauss].normal[0] = interface.normal[0];
				flux[istage].gauss[igauss].normal[1] = interface.normal[1];
				flux[istage].gauss[igauss].length = interface.length;
			}
		}
	}
}
void Set_Gauss_Coordinate(Interface2d& xinterface, Interface2d& yinterface)
{
	for (int num_guass = 0; num_guass < gausspoint; num_guass++)
	{
		// first is gauss parameter
		xinterface.gauss[num_guass].x = xinterface.x;
		xinterface.gauss[num_guass].y = xinterface.y + gauss_loc[num_guass] * 0.5 * xinterface.length;
		xinterface.gauss[num_guass].normal[0] = xinterface.normal[0];
		xinterface.gauss[num_guass].normal[1] = xinterface.normal[1];
		// each gauss point contain left center right point
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].left, num_guass, xinterface.length);
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].right, num_guass, xinterface.length);
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].center, num_guass, xinterface.length);

		yinterface.gauss[num_guass].x = yinterface.x + gauss_loc[num_guass] * 0.5 * yinterface.length;
		yinterface.gauss[num_guass].y = yinterface.y;
		yinterface.gauss[num_guass].normal[0] = yinterface.normal[0];
		yinterface.gauss[num_guass].normal[1] = yinterface.normal[1];
		// each gauss point contain left center right point
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].left, num_guass, yinterface.length);
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].right, num_guass, yinterface.length);
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].center, num_guass, yinterface.length);

	}
}
void Set_Gauss_Intergation_Location_x(Point2d& xgauss, int index, double h)
{
	xgauss.x = gauss_loc[index] * h / 2.0;
}
void Set_Gauss_Intergation_Location_y(Point2d& ygauss, int index, double h)
{
	ygauss.x = gauss_loc[index] * h / 2.0;
}


void Update(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage, Interface2d* xinterfaces, Interface2d* yinterfaces)
{

	if (stage > block.stages)
	{
		cout << "wrong middle stage,pls check the time marching setting" << endl;
		exit(0);
	}

	double dt = block.dt;

#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				for (int var = 0; var < 4; var++)
				{
					double Flux = 0.0;
					for (int k = 0; k < stage + 1; k++)
					{

						Flux = Flux
							+ gauss_weight[num_gauss] *
							(block.timecoefficient[stage][k][0] * xfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].f[var]
							+ block.timecoefficient[stage][k][1] * xfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].derf[var]);

					}
					xfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] = Flux * dt;
					// calculate the final flux of the interface, in x in xfluxes, by the obtained the flux and its derivative (f, def, der2f) at guass points, in xfluxes, and the corresponding weight factors
					// calculate by several stages according to the time marching method. same for yfluxes
				}
			}
		}
	}

#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost + 1; j++)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				for (int var = 0; var < 4; var++)
				{
					double Flux = 0.0;
					for (int k = 0; k < stage + 1; k++)
					{
						Flux = Flux
							+ gauss_weight[num_gauss] *
							(block.timecoefficient[stage][k][0] * yfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].f[var]
							+ block.timecoefficient[stage][k][1] * yfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].derf[var]);
					}
					yfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] = Flux * dt;
				}
			}
		}
	}

	
	Update_with_gauss(fluids, xfluxes, yfluxes, block, stage);
}
void Update_with_gauss(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage)
{
	//Note : calculate the final flux of the cell, in fluids, by the obtained final flux of interface, in xfluxes and yfluxes
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost + 1; j++)
		{
			int face = i * (block.ny + 1) + j;
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				Local_to_Global(yfluxes[face][stage].gauss[num_gauss].x, yfluxes[face][stage].gauss[num_gauss].normal);
				Local_to_Global(xfluxes[face][stage].gauss[num_gauss].x, xfluxes[face][stage].gauss[num_gauss].normal);
			}
		}
	}

//#pragma omp parallel  for
	//for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	//{
	//	for (int j = block.ghost; j < block.nodey + block.ghost; j++)
	//	{
	//		int cell = i * (block.ny) + j;
	//		int face = i * (block.ny + 1) + j;

	//		for (int var = 0; var < 4; var++)
	//		{
	//			fluids[cell].convar[var] = fluids[cell].convar_old[var]; //get the Wn from convar_old
	//			double total_flux = 0.0;
	//			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	//			{
	//				total_flux += yfluxes[face][stage].gauss[num_gauss].length * yfluxes[face][stage].gauss[num_gauss].x[var];
	//				total_flux += -yfluxes[face + 1][stage].gauss[num_gauss].length * yfluxes[face + 1][stage].gauss[num_gauss].x[var];
	//				total_flux += xfluxes[face][stage].gauss[num_gauss].length * xfluxes[face][stage].gauss[num_gauss].x[var];
	//				total_flux += -xfluxes[face + block.ny + 1][stage].gauss[num_gauss].length * xfluxes[face + block.ny + 1][stage].gauss[num_gauss].x[var];
	//			}
	//			fluids[cell].convar[var] += total_flux / fluids[cell].area;

	//		}
	//	}

	//}

	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				for (int var = 0; var < 4; var++)
				{
					double local;
					int cell = i * (block.ny) + j;
					int face = i * (block.ny + 1) + j;

					//fluids[cell].convar[var] = fluids[cell].convar_old[var]; //get the Wn from convar_old

					fluids[i * block.ny + j].dfdx[var] = -(xfluxes[(i + 1) * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] - xfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var]) / block.dx; // y0 du/dx
					fluids[i * block.ny + j].dfdy[var] = -(yfluxes[i * (block.ny + 1) + j + 1][stage].gauss[num_gauss].x[var] - yfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var]) / block.dy;   // y0 du/dy
					fluids[i * block.ny + j].df[var] = fluids[i * block.ny + j].dfdx[var] + fluids[i * block.ny + j].dfdy[var];   // y0  du/dx + du/dy

					//fluids[i * block.ny + j].convar[var] += fluids[i * block.ny + j].df[var];
					fluids[i * block.ny + j].convar[var] = fluids[cell].convar_old[var] + fluids[i * block.ny + j].df[var];

				}
				
			}
		}
	}
}


void Check_Order_Reduce_by_Lambda_2D(bool& order_reduce, double* convar)
{
	order_reduce = false;
	double lambda;
	lambda = Lambda(convar[0], convar[1] / convar[0], convar[2] / convar[0], convar[3]);
	//if lambda <0, then reduce to the first order
	if (lambda <= 0.0 || (lambda == lambda) == false)
	{
		order_reduce = true;
	}
}
void Recompute_KFVS_1st(Fluid2d& fluid, double* center, double* left, double* right, double* down, double* up, Block2d& block)
{
	double flux_left[4];
	double flux_right[4];
	double flux_down[4];
	double flux_up[4];
	KFVS_1st(flux_left, left, center, block.dt);
	KFVS_1st(flux_right, center, right, block.dt);

	double center_tmp[4], up_tmp[4], down_tmp[4];
	YchangetoX(center_tmp, center);
	YchangetoX(up_tmp, up);
	YchangetoX(down_tmp, down);

	KFVS_1st(flux_down, down_tmp, center_tmp, block.dt);
	KFVS_1st(flux_up, center_tmp, up_tmp, block.dt);

	fluid.convar[0] = fluid.convar_old[0] + 1.0 / block.dy * (flux_down[0] - flux_up[0]);
	fluid.convar[1] = fluid.convar_old[1] + 1.0 / block.dy * (-flux_down[2] + flux_up[2]);
	fluid.convar[2] = fluid.convar_old[2] + 1.0 / block.dy * (flux_down[1] - flux_up[1]);
	fluid.convar[3] = fluid.convar_old[3] + 1.0 / block.dy * (flux_down[3] - flux_up[3]);

	fluid.convar[0] = fluid.convar[0] + 1.0 / block.dx * (flux_left[0] - flux_right[0]);
	fluid.convar[1] = fluid.convar[1] + 1.0 / block.dx * (flux_left[1] - flux_right[1]);
	fluid.convar[2] = fluid.convar[2] + 1.0 / block.dx * (flux_left[2] - flux_right[2]);
	fluid.convar[3] = fluid.convar[3] + 1.0 / block.dx * (flux_left[3] - flux_right[3]);
}
void KFVS_1st(double* flux, double* left, double* right, double dt)
{
	double density_left, density_right;

	double u_left, u_right, v_left, v_right;
	//λ
	double lambda_left, lambda_right;

	density_left = left[0];
	density_right = right[0];

	u_left = U(density_left, left[1]);
	v_left = V(density_left, left[2]);
	lambda_left = Lambda(density_left, u_left, v_left, left[3]);

	u_right = U(density_right, right[1]);
	v_right = V(density_right, right[2]);
	lambda_right = Lambda(density_right, u_right, v_right, right[3]);

	MMDF1st m2(u_left, v_left, lambda_left);
	MMDF1st m3(u_right, v_right, lambda_right);

	//t4part
	flux[0] =
		density_left * dt * m2.uplus[1] +
		density_right * dt * m3.uminus[1];

	flux[1] =
		density_left * dt * m2.uplus[2] +
		density_right * dt * m3.uminus[2];

	flux[2] =
		density_left * dt * m2.uplus[1] * m2.vwhole[1] +
		density_right * dt * m3.uminus[1] * m3.vwhole[1];

	flux[3] =
		density_left * 0.5 * dt * (m2.uplus[3] + m2.uplus[1] * m2.vwhole[2] + m2.uplus[1] * m2.xi2) +
		density_right * 0.5 * dt * (m3.uminus[3] + m3.uminus[1] * m3.vwhole[2] + m3.uminus[1] * m3.xi2);
}



void SetNonUniformMesh(Block2d& block, Fluid2d* fluids,
	Interface2d* xinterfaces, Interface2d* yinterfaces,
	Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, double** node2d)
{
	block.xcell_begin = block.ghost;
	block.xcell_end = block.ghost + block.nodex - 1;
	block.ycell_begin = block.ghost;
	block.ycell_end = block.ghost + block.nodey - 1;

	block.xinterface_begin_n = block.ghost;
	block.xinterface_end_n = block.ghost + block.nodex;
	block.xinterface_begin_t = block.ghost;
	block.xinterface_end_t = block.ghost + block.nodex - 1;

	block.yinterface_begin_n = block.ghost;
	block.yinterface_end_n = block.ghost + block.nodey;
	block.yinterface_begin_t = block.ghost;
	block.yinterface_end_t = block.ghost + block.nodey - 1;

#pragma omp parallel  for
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			int cindex = i * block.ny + j;
			int node_index[4];
			node_index[0] = i * (block.ny + 1) + j;
			node_index[1] = (i + 1) * (block.ny + 1) + j;
			node_index[2] = (i + 1) * (block.ny + 1) + j + 1;
			node_index[3] = (i) * (block.ny + 1) + j + 1;
			Set_cell_geo_from_quad_node(fluids[cindex],
				node2d[node_index[0]], node2d[node_index[1]],
				node2d[node_index[2]], node2d[node_index[3]]);
		}
	}

#pragma omp parallel  for
	for (int i = 0; i < block.nx + 1; i++)
	{
		for (int j = 0; j < block.ny + 1; j++)
		{
			int interface_index = i * (block.ny + 1) + j;
			int node_index[3];
			node_index[0] = i * (block.ny + 1) + j;
			node_index[1] = (i + 1) * (block.ny + 1) + j;
			node_index[2] = (i) * (block.ny + 1) + j + 1;
			if (j < block.ny)
			{
				Set_interface_geo_from_two_node
				(xinterfaces[interface_index], node2d[node_index[0]], node2d[node_index[2]], 0);
			}
			if (i < block.nx)
			{
				Set_interface_geo_from_two_node
				(yinterfaces[interface_index], node2d[node_index[0]], node2d[node_index[1]], 1);
			}

			Copy_geo_from_interface_to_line(xinterfaces[interface_index]);
			Copy_geo_from_interface_to_line(yinterfaces[interface_index]);

			xinterfaces[interface_index].gauss = new Recon2d[gausspoint];
			yinterfaces[interface_index].gauss = new Recon2d[gausspoint];
			if (j < block.ny)
			{
				Set_Gauss_Coordinate_general_mesh_x(xinterfaces[interface_index],
					node2d[node_index[0]], node2d[node_index[2]]);
			}
			if (i < block.nx)
			{
				Set_Gauss_Coordinate_general_mesh_y(yinterfaces[interface_index],
					node2d[node_index[0]], node2d[node_index[1]]);
			}


			Copy_geo_from_interface_to_flux
			(xinterfaces[interface_index], xfluxes[interface_index], block.stages);
			Copy_geo_from_interface_to_flux
			(yinterfaces[interface_index], yfluxes[interface_index], block.stages);

		}
	}

}


void Set_cell_geo_from_quad_node
(Fluid2d& fluid, double* n1, double* n2, double* n3, double* n4)
{
	//4            //3
	/////////////////
	//             //
	//             //
	//             //
	//             //
	//             //
	/////////////////
	//1            //2
	Copy_Array(&fluid.node[0], n1, 2);
	Copy_Array(&fluid.node[2], n2, 2);
	Copy_Array(&fluid.node[4], n3, 2);
	Copy_Array(&fluid.node[6], n4, 2);
	fluid.coordx = (n1[0] + n2[0] + n4[0] + n3[0]) / 4.0;
	fluid.coordy = (n1[1] + n2[1] + n4[1] + n3[1]) / 4.0;
	fluid.area = 0.5 * sqrt((pow((n2[0] - n1[0]), 2) + pow((n2[1] - n1[1]), 2)) *
		(pow((n4[0] - n1[0]), 2) + pow((n4[1] - n1[1]), 2)) -
		pow(((n2[0] - n1[0]) * (n4[0] - n1[0]) +
			(n2[1] - n1[1]) * (n4[1] - n1[1])), 2)) +
		0.5 * sqrt((pow((n2[0] - n3[0]), 2) + pow((n2[1] - n3[1]), 2)) *
			(pow((n4[0] - n3[0]), 2) + pow((n4[1] - n3[1]), 2)) -
			pow(((n2[0] - n3[0]) * (n4[0] - n3[0]) +
				(n2[1] - n3[1]) * (n4[1] - n3[1])), 2));
	fluid.dx = sqrt(pow((0.5 * (n2[0] + n3[0]) - 0.5 * (n1[0] + n4[0])), 2)
		+ pow((0.5 * (n3[1] - n4[1]) + 0.5 * (n2[1] - n1[1])), 2));
	fluid.dy = sqrt(pow((0.5 * (n4[0] + n3[0]) - 0.5 * (n1[0] + n2[0])), 2) +
		pow((0.5 * (n4[1] + n3[1]) - 0.5 * (n1[1] + n2[1])), 2));

}


void Set_interface_geo_from_two_node
(Interface2d& interface, double* n1, double* n2, int direction)
{
	interface.x = (n1[0] + n2[0]) / 2.0;
	interface.y = (n1[1] + n2[1]) / 2.0;
	interface.length =
		sqrt(pow((n2[0] - n1[0]), 2) + pow((n2[1] - n1[1]), 2));
	if (direction == 0) //we see it as x direction
	{
		interface.normal[0] =
			(n2[1] - n1[1]) / interface.length;
		interface.normal[1] =
			-(n2[0] - n1[0]) / interface.length;
	}
	else
	{
		interface.normal[0] =
			-(n2[1] - n1[1]) / interface.length;
		interface.normal[1] =
			(n2[0] - n1[0]) / interface.length;
	}

}


void Set_Gauss_Coordinate_general_mesh_x
(Interface2d& xinterface, double* node0, double* nodex1)
{
	double xoff[2];
	xoff[0] = nodex1[0] - node0[0]; xoff[1] = nodex1[1] - node0[1];
	for (int num_guass = 0; num_guass < gausspoint; num_guass++)
	{
		// first is gauss parameter
		xinterface.gauss[num_guass].x = xinterface.x + gauss_loc[num_guass] * 0.5 * xoff[0];
		xinterface.gauss[num_guass].y = xinterface.y + gauss_loc[num_guass] * 0.5 * xoff[1];
		xinterface.gauss[num_guass].normal[0] = xinterface.normal[0];
		xinterface.gauss[num_guass].normal[1] = xinterface.normal[1];
		// each gauss point contain left center right point 
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].left, num_guass, xinterface.length);
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].right, num_guass, xinterface.length);
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].center, num_guass, xinterface.length);
	}
}

void Set_Gauss_Coordinate_general_mesh_y
(Interface2d& yinterface, double* node0, double* nodey1)
{
	double yoff[2];
	yoff[0] = nodey1[0] - node0[0]; yoff[1] = nodey1[1] - node0[1];
	for (int num_guass = 0; num_guass < gausspoint; num_guass++)
	{
		yinterface.gauss[num_guass].x = yinterface.x + gauss_loc[num_guass] * 0.5 * yoff[0];
		yinterface.gauss[num_guass].y = yinterface.y + gauss_loc[num_guass] * 0.5 * yoff[1];
		yinterface.gauss[num_guass].normal[0] = yinterface.normal[0];
		yinterface.gauss[num_guass].normal[1] = yinterface.normal[1];
		// each gauss point contain left center right point 
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].left, num_guass, yinterface.length);
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].right, num_guass, yinterface.length);
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].center, num_guass, yinterface.length);
	}
}


void Residual2d(Fluid2d* fluids, Block2d block, int outputstep)
{
	if (block.step % outputstep == 0)
	{
		int order = block.ghost;
		double residual[4];
		double sum_old[4];
		for (int k = 0; k < 4; k++)
		{
			residual[k] = 0.0;
			sum_old[k] = 0.0;
		}
		for (int i = order; i < block.nx - order; i++)
		{
			for (int j = order; j < block.ny - order; j++)
			{
				int index = i * block.ny + j;
				for (int k = 0; k < 4; k++)
				{
					residual[k] = residual[k] + abs(fluids[index].convar[k] - fluids[index].convar_old[k]);
					sum_old[k] = sum_old[k] + abs(fluids[index].convar_old[k]);
				}
			}
		}
		for (int k = 0; k < 4; k++)
		{
			residual[k] = residual[k] / (sum_old[k] + 1e-15) / block.dt;
		}

		cout << "step= " << block.step
			<< " density " << residual[0] << " u " << residual[1]
			<< " v " << residual[2] << " densityE " << residual[3] << endl;


		ofstream out;
		if (block.step == outputstep)
		{
			out.open("result/residual-2D.plt", ios::out);
		}
		else
		{
			out.open("result/residual-2D.plt", ios::ate | ios::out | ios::in);
		}

		if (!out.is_open())
		{
			cout << "cannot find residual-2D.plt" << endl;
			cout << "a new case will start" << endl;
		}

		if (block.step == outputstep)
		{
			out << "# CFL number is " << block.CFL << endl;
			out << "# tau type (0 refer euler, 1 refer NS) is " << tau_type << endl;
			out << "# time marching strategy is " << block.stages << " stage method" << endl;
			out << "variables = step,density_residual,m_residual,n_residual,E_residual" << endl;

		}

		//output the data
		out << block.step << " "
			<< log10(residual[0]) << " " << log10(residual[1]) << " "
			<< log10(residual[2]) << " " << log10(residual[3]) << endl;
		out.close();
	}
}
