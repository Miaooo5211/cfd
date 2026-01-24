#include"cylinder.h"
#include"boundary_layer.h"
void output2d(Fluid2d* fluids, Block2d block);
void cylinder()
{
	Runtime runtime;
	runtime.start_initial = clock();
	Block2d block;
	block.uniform = false;

	block.ghost = 9;

	double tstop = 0.2; //Sod
	block.CFL = 1.0;

	K = 3;
	Gamma = 1.4;
	Pr = 1.0;

	double renum = 1e5;
	double den_ref = 1.0;
	double u_ref = 0.15;
	double L = 100;

	//prepare the boundary condtion function
	Fluid2d* bcvalue = new Fluid2d[1];
	bcvalue[0].primvar[0] = den_ref;
	bcvalue[0].primvar[1] = u_ref;
	bcvalue[0].primvar[2] = 0.0;
	bcvalue[0].primvar[3] = den_ref / Gamma;

	Mu = den_ref * u_ref * L / renum;
	cout << Mu << endl;

	gausspoint = 1;
	SetGuassPoint();
	reconstruction_variable = characteristic;
	wenotype = teno;

	cellreconstruction_2D_normal = TENO7_normal;


	flux_function_2d = CGKS2D;

	timecoe_list_2d = RK2_2D;
	Initial_stages(block);
	block.stages = 5;


	Fluid2d* fluids = NULL;
	Interface2d* xinterfaces = NULL;
	Interface2d* yinterfaces = NULL;
	Flux2d_gauss** xfluxes = NULL;
	Flux2d_gauss** yfluxes = NULL;

	string grid = add_mesh_directory_modify_for_linux()
		+ "structured-mesh/blast200.plt";

	Read_structured_mesh
	(grid, &fluids, &xinterfaces, &yinterfaces, &xfluxes, &yfluxes, block);

	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			int index = i * block.ny + j;

			if (fluids[index].coordx < 0.5)
			{
				fluids[index].convar[0] = 1.0;
				fluids[index].convar[1] = 0.0;
				fluids[index].convar[2] = 0.0;
				fluids[index].convar[3] = 2.5;
			}
			else
			{
				fluids[index].convar[0] = 0.125;
				fluids[index].convar[1] = 0.0;
				fluids[index].convar[2] = 0.0;
				fluids[index].convar[3] = 0.25;
			}

			for (int var = 0; var < 4; var++)
			{
				fluids[index].D_u[var] = fluids[index].convar[var];

				fluids[index].convar0[var] = fluids[index].convar[var];
				fluids[index].convar1[var] = fluids[index].convar[var];
				fluids[index].convar2[var] = fluids[index].convar[var];
				fluids[index].convar3[var] = fluids[index].convar[var];
				fluids[index].convar4[var] = fluids[index].convar[var];

			}
		}
	}



	runtime.finish_initial = clock();
	block.t = 0;
	block.step = 0;
	int inputstep = 1;

	block.dx = 1.0 / block.nodex;
	block.dy = 1.0 / block.nodey;

	while (block.t < tstop)
	{
		if (block.step % inputstep == 0)
		{
			cout << "pls cin interation step, if input is 0, then the program will exit " << endl;
			cin >> inputstep;
			if (inputstep == 0)
			{
				output2d(fluids, block);
				break;
			}
		}
		if (runtime.start_compute == 0.0)
		{
			runtime.start_compute = clock();
			cout << "runtime-start " << endl;
		}


		CopyFluid_new_to_old(fluids, block);
		Convar_to_Primvar(fluids, block);


		block.dt = Get_CFL1(block, fluids, tstop);

		//-------------------DeC7_A--------------------
		int iter = 7;
		cellreconstruction_2D_normal = TENO7_normal;

		for (int k = 0; k < iter; k++)  //correction
		{
			if (k == 0)
			{
				boundaryforBoundary_layerSod(fluids, block, bcvalue[0]);

				Reconstruction_within_cell0(xinterfaces, yinterfaces, fluids, block, 0);

				Calculate_flux(xfluxes, yfluxes, xinterfaces, yinterfaces, block, 0, fluids);

				comput_du_dt(fluids, xfluxes, yfluxes, block, 0, xinterfaces, yinterfaces);

				Update_DeC7_1(fluids, xfluxes, yfluxes, block, xinterfaces, yinterfaces);
			}
			else
			{
				for (int i = 1; i < block.stages; i++) //midstep
				{
					convarset(fluids, block, bcvalue[0], i);

					boundaryforBoundary_layerSod(fluids, block, bcvalue[0]);

					Reconstruction_within_cell1(xinterfaces, yinterfaces, fluids, block);

					Calculate_flux(xfluxes, yfluxes, xinterfaces, yinterfaces, block, i, fluids);

					comput_du_dt(fluids, xfluxes, yfluxes, block, i, xinterfaces, yinterfaces);

					Duset(fluids, block, bcvalue[0], i);
				}
				Update_DeC7(fluids, xfluxes, yfluxes, block, xinterfaces, yinterfaces);
			}
		}
		finalset(fluids, block, bcvalue[0]);

		//-------------------DeC7--------------------
		//int iter = 7;
		//cellreconstruction_2D_normal = TENO70_normal;

		//for (int k = 0; k < iter; k++)  //correction
		//{
		//	if (k == 0)
		//	{
		//		boundaryforBoundary_layerSod(fluids, block, bcvalue[0]);

		//		Reconstruction_within_cell(xinterfaces, yinterfaces, fluids, block);

		//		Calculate_flux(xfluxes, yfluxes, xinterfaces, yinterfaces, block, 0, fluids);

		//		comput_du_dt(fluids, xfluxes, yfluxes, block, 0, xinterfaces, yinterfaces);

		//		Update_DeC7_1(fluids, xfluxes, yfluxes, block, xinterfaces, yinterfaces);
		//	}
		//	else
		//	{
		//		for (int i = 1; i < block.stages; i++) //midstep
		//		{
		//			convarset(fluids, block, bcvalue[0], i);

		//			boundaryforBoundary_layerSod(fluids, block, bcvalue[0]);

		//			Reconstruction_within_cell(xinterfaces, yinterfaces, fluids, block);

		//			Calculate_flux(xfluxes, yfluxes, xinterfaces, yinterfaces, block, i, fluids);

		//			comput_du_dt(fluids, xfluxes, yfluxes, block, i, xinterfaces, yinterfaces);

		//			Duset(fluids, block, bcvalue[0], i);
		//		}
		//		Update_DeC7(fluids, xfluxes, yfluxes, block, xinterfaces, yinterfaces);
		//	}
		//}
		//finalset(fluids, block, bcvalue[0]);

		////---------------- RK44----------------
		//block.stages = 4;
		//cellreconstruction_2D_normal = TENO70_normal;

		//for (int i = 0; i < block.stages; i++)
		//{
		//	boundaryforBoundary_layerSod(fluids, block, bcvalue[0]);

		//	Reconstruction_within_cell(xinterfaces, yinterfaces, fluids, block);

		//	Calculate_flux(xfluxes, yfluxes, xinterfaces, yinterfaces, block, i, fluids);

		//	Update_RK44(fluids, xfluxes, yfluxes, block, i, xinterfaces, yinterfaces);
		//}

		++block.step;
		block.t = block.t + block.dt;

		if (block.step % 10 == 0)
		{
			cout << "step 10 time is " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;

		}
		Residual2d(fluids, block, 10);


	}
	runtime.finish_compute = clock();
	;
	cout << "\n the total running time is " << (double)(runtime.finish_compute - runtime.start_initial) / CLOCKS_PER_SEC << "Ãë£¡" << endl;
	cout << "\n the time for initializing is " << (double)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << "Ãë£¡" << endl;
	cout << "\n the time for computing is " << (double)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << "Ãë£¡" << endl;

	output2d(fluids, block);

	int order = block.ghost;
	int i;


	stringstream name;

	name << "result/line-dentCFL=" << block.CFL << ".plt" << endl;


	string s;
	name >> s;
	ofstream out(s);
	out << "variables =x,density,us,vs,pressure" << endl;
	out << "zone i = " << 1 << ",j = " << block.nodex << ", F=POINT" << endl;


	for (int j = order; j < block.nx - order; j++)
	{
		int index = j * block.ny + order;
		Convar_to_primvar_2D(fluids[index].primvar, fluids[index].convar);

		out << fluids[index].coordx << " "
			<< fluids[index].primvar[0] << " "
			<< fluids[index].primvar[1] << " "
			<< fluids[index].primvar[2] << " "
			<< fluids[index].sensor << " "
			<< endl;
	}

	out.close();


}


void finalset(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
#pragma omp parallel for collapse(2)
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			int index = i * block.ny + j;

			for (int var = 0; var < 4; var++)
			{
				fluids[index].convar_old[var] = fluids[index].convar0[var];
				fluids[index].convar[var] = fluids[index].convar4[var];

				fluids[index].convar0[var] = fluids[index].convar[var];
				fluids[index].convar1[var] = fluids[index].convar[var];
				fluids[index].convar2[var] = fluids[index].convar[var];
				fluids[index].convar3[var] = fluids[index].convar[var];
				fluids[index].convar4[var] = fluids[index].convar[var];
			}
		}
	}
}

void convarset(Fluid2d* fluids, Block2d block, Fluid2d bcvalue, int stage)
{
#pragma omp parallel for collapse(2)
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			int index = i * block.ny + j;

			for (int var = 0; var < 4; var++)
			{
				if (stage == 0)
				{
					fluids[index].convar[var] = fluids[index].convar0[var];
				}
				else if (stage == 1)
				{
					fluids[index].convar[var] = fluids[index].convar1[var];
				}
				else if (stage == 2)
				{
					fluids[index].convar[var] = fluids[index].convar2[var];
				}
				else if (stage == 3)
				{
					fluids[index].convar[var] = fluids[index].convar3[var];
				}
				else
				{
					fluids[index].convar[var] = fluids[index].convar4[var];
				}
			}
		}
	}
}

void Duset(Fluid2d* fluids, Block2d block, Fluid2d bcvalue, int stage)
{
#pragma omp parallel for collapse(2)
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			int index = i * block.ny + j;

			for (int var = 0; var < 4; var++)
			{
				if (stage == 0)
				{
					fluids[index].D_du_ut0[var] = fluids[index].D_du_ut[var];
				}
				else if (stage == 1)
				{
					fluids[index].D_du_ut1[var] = fluids[index].D_du_ut[var];
				}
				else if (stage == 2)
				{
					fluids[index].D_du_ut2[var] = fluids[index].D_du_ut[var];
				}
				else if (stage == 3)
				{
					fluids[index].D_du_ut3[var] = fluids[index].D_du_ut[var];
				}
				else
				{
					fluids[index].D_du_ut4[var] = fluids[index].D_du_ut[var];
				}
			}
		}
	}
}

void boundaryforBoundary_layerSod(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{

	for (int i = block.ghost - 1; i >= 0; i--)
	{
		for (int j = 0; j < block.ny; j++)
		{
			fluids[i * block.ny + j].convar[0] = 1.0;
			fluids[i * block.ny + j].convar[1] = 0.0;
			fluids[i * block.ny + j].convar[2] = 0.0;
			fluids[i * block.ny + j].convar[3] = 2.5;
		}
	}

	//outlet-right

	for (int i = block.nx - block.ghost; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			fluids[i * block.ny + j].convar[0] = 0.125;
			fluids[i * block.ny + j].convar[1] = 0.0;
			fluids[i * block.ny + j].convar[2] = 0.0;
			fluids[i * block.ny + j].convar[3] = 0.25;
		}
	}

	//free-down

	for (int j = block.ghost - 1; j >= 0; j--)
	{
		for (int i = 0; i < block.nx; i++)
		{
			fluids[i * block.ny + j].convar[0] = fluids[i * block.ny + j + 1].convar[0];
			fluids[i * block.ny + j].convar[1] = fluids[i * block.ny + j + 1].convar[1];
			fluids[i * block.ny + j].convar[2] = fluids[i * block.ny + j + 1].convar[2];
			fluids[i * block.ny + j].convar[3] = fluids[i * block.ny + j + 1].convar[3];
		}
	}


	//free-up

	for (int j = block.ny - block.ghost; j < block.ny; j++)
	{
		for (int i = 0; i < block.nx; i++)
		{
			fluids[i * block.ny + j].convar[0] = fluids[i * block.ny + j - 1].convar[0];
			fluids[i * block.ny + j].convar[1] = fluids[i * block.ny + j - 1].convar[1];
			fluids[i * block.ny + j].convar[2] = fluids[i * block.ny + j - 1].convar[2];
			fluids[i * block.ny + j].convar[3] = fluids[i * block.ny + j - 1].convar[3];
		}
	}
}



