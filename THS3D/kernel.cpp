#include "kernel.h"
#include <ctime>
 #include <numeric>

// void kernel::run()
// {
// 	int niter(20);
// 	unsigned int i;
// 	TruncatedTspline_3D tt3;
// 	//tt3.SetProblem("../io/hex_input/cube_dense");
// 	//tt3.SetProblem("../io/hex_input/cube4");
// 	tt3.SetProblem("../io/hex_input/cube_coarse_0");
// 	string fld("../io/hex_out_2/"), fn("test16_");
// 	vector<int> dof_list(niter,0);
// 	vector<double> err_list(niter,0.);
// 	for (int itr = 0; itr < niter; itr++)
// 	{
// 		vector<BezierElement3D> bzmesh;
// 		vector<int> IDBC;
// 		vector<double> gh, err;
// 		//cout << "interface...\n";
// 		tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
// 		//tt3.AnalysisInterface_Laplace(bzmesh, IDBC, gh);
// 		//cout << "interface done!\n";
// 		//getchar();

// 		//cout << "npt: " << IDBC.size() << "\n";
// 		//getchar();

// 		Laplace lap;
// 		lap.SetProblem(IDBC, gh);
// 		stringstream ss;
// 		ss << itr;
// 		//cout << "before simulation...\n";
// 		//getchar();
// 		lap.Run(bzmesh, fld+fn+ss.str(), err);
// 		//cout << "simulation done!\n";
// 		//getchar();

// 		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

// 		double errL2(0.);
// 		for (i = 0; i < err.size(); i++) errL2 += err[i];
// 		errL2 = sqrt(errL2);
// 		dof_list[itr] = IDBC.size();
// 		err_list[itr] = errL2;
// 		//cout << "DOF: " << IDBC.size() << "\n";
// 		//cout << "L2-norm error: " << errL2 << "\n";
// 		//getchar();

// 		if (itr < niter - 1)
// 		{
// 			cout << itr << " refining...\n";
// 			//distribute error
// 			vector<array<double, 2>> eh(bzmesh.size());
// 			for (i = 0; i < bzmesh.size(); i++)
// 			{
// 				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
// 			}
// 			vector<array<int, 2>> rfid, gst;
// 			//tt3.Identify_Poisson(eh, err, rfid, gst);
// 			tt3.Identify_Poisson_1(eh, err, rfid, gst);
// 			//tt3.Identify_Laplace(eh, err, rfid, gst);
// 			tt3.Refine(rfid, gst);
// 			cout << "Refining done\n";
// 			tt3.OutputGeom_All(fld+fn + ss.str()+"_geom");
// 			//cout << "Output Geom done!\n";
// 			//getchar();
// 		}
// 	}

// 	//output error
// 	output_err(fld + fn + "err", dof_list, err_list);
// }

// void kernel::run_complex()
// {
// 	int niter(2);
// 	unsigned int i;
// 	TruncatedTspline_3D tt3;

// 	tt3.SetProblem("../io/hex_input/statue");

// 	string fld("../io/complex2/");

// 	string fn("statue1_");

// 	double remove[3][2];
// 	tt3.GetRemoveRegion(remove);

// 	vector<int> dof_list(niter, 0);
// 	vector<double> err_list(niter, 0.);
// 	vector<array<int, 2>> ebc, pbc;
// 	vector<double> edisp, pdisp;
// 	tt3.SetInitialBC(ebc, edisp, pbc, pdisp);
// 	for (int itr = 0; itr < niter; itr++)
// 	{
// 		vector<BezierElement3D> bzmesh;
// 		vector<int> IDBC;
// 		vector<double> gh, err;
// 		tt3.SetBC(ebc, edisp, pbc, pdisp);
// 		//cout << pbc.size() << "\n";
// 		//cout << "interface...\n";
// 		//tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
// 		//tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
// 		tt3.AnalysisInterface_Laplace(pbc,pdisp,bzmesh, IDBC, gh);
// 		//cout << "interface done!\n";
// 		//getchar();

// 		//cout << "npt: " << IDBC.size() << "\n";
// 		//getchar();

// 		Laplace lap;
// 		lap.SetProblem(IDBC, gh);
// 		stringstream ss;
// 		ss << itr;
// 		//cout << "before simulation...\n";
// 		//getchar();
// 		lap.GetRemoveRegion(remove);
// 		lap.Run(bzmesh, fld + fn + ss.str(), err);
// 		//cout << "simulation done!\n";
// 		//getchar();

// 		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

// 		double errL2(0.);
// 		for (i = 0; i < err.size(); i++) errL2 += err[i];
// 		errL2 = sqrt(errL2);
// 		dof_list[itr] = IDBC.size();
// 		err_list[itr] = errL2;
// 		//cout << "DOF: " << IDBC.size() << "\n";
// 		//cout << "L2-norm error: " << errL2 << "\n";
// 		//getchar();

// 		output_err(fld + fn+ss.str() + "_err", dof_list, err_list);

// 		if (itr < niter - 1)
// 		{
// 			cout << itr << " refining...\n";
// 			//distribute error
// 			vector<array<double, 2>> eh(bzmesh.size());
// 			for (i = 0; i < bzmesh.size(); i++)
// 			{
// 				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
// 			}
// 			vector<array<int, 2>> rfid, gst;
// 			//tt3.Identify_Poisson_1(eh, err, rfid, gst);
// 			tt3.Identify_Laplace(eh, err, rfid, gst);
// 			//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
// 			//tt3.InputRefineID("../io/benchmark1/cube1_" + ss.str(), rfid, gst);
// 			tt3.Refine(rfid, gst);
// 			cout << "Refining done\n";
// 			//tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
// 			//cout << "Output Geom done!\n";
// 			//getchar();
// 		}
// 	}

// 	//output error
// 	//output_err(fld + fn + "err", dof_list, err_list);
// }


void kernel::OutputMesh(const vector<BezierElement3D>& bzmesh, string fn)
{
	int cn[8] = { 0, 3, 15, 12, 48, 51, 63, 60 };
	string fname = fn + "bzmesh.vtk";
	ofstream fout;
	fout.open(fname.c_str());	

	// cout << "ck2" << endl;

	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nBezier mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << 8 * bzmesh.size() << " float\n";		
		for (int i = 0; i<bzmesh.size(); i++)
		{
			for (int j = 0; j < 8; j++)
			{
				fout << bzmesh[i].pts[cn[j]][0] << " " << bzmesh[i].pts[cn[j]][1] << " " << bzmesh[i].pts[cn[j]][2] << "\n";
			}
		}
		fout << "\nCELLS " << bzmesh.size() << " " << 9 * bzmesh.size() << '\n';
		for (int i = 0; i<bzmesh.size(); i++)
		{
			fout << "8 " << 8 * i << " " << 8 * i + 1 << " " << 8 * i + 2 << " " << 8 * i + 3
				<< " " << 8 * i + 4 << " " << 8 * i + 5 << " " << 8 * i + 6 << " " << 8 * i + 7 << '\n';
		}
		fout << "\nCELL_TYPES " << bzmesh.size() << '\n';
		for (int i = 0; i<bzmesh.size(); i++)
		{
			fout << "12\n";
		}
		//fout << "POINT_DATA " << sdisp.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		//for (int i = 0; i<sdisp.size(); i++)
		//{
		//	fout << sdisp[i] << "\n";
		//}
		fout << "\nCELL_DATA " << bzmesh.size() << "\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		for (int i = 0; i < bzmesh.size(); i++)
		{
			fout << bzmesh[i].type << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
	
	// cout << "ck3" << endl;

	string fname3(fn + "bzmeshinfo.txt");
	//ofstream fout;
	fout.open(fname3.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() << "\n";
		for (int i = 0; i<bzmesh.size(); i++)
		{
			for (int l = 0; l < bzmesh[i].IEN.size(); l++)
			{
				fout << bzmesh[i].IEN[l] + 1;
				if (l == bzmesh[i].IEN.size() - 1)
				{
					fout << "\n";
				}
				else
				{
					fout << " ";
				}
			}
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname3 << '\n';
	}

	// cout << "ck4" << endl;

	string fname1(fn + "cmat.txt");
	//ofstream fout;
	fout.open(fname1.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() << "\n";
		for (int i = 0; i<bzmesh.size(); i++)
		{
			fout << i << " " << bzmesh[i].IEN.size() << " " << bzmesh[i].type << "\n";
			for (int l = 0; l < bzmesh[i].IEN.size(); l++)
			{
				fout << bzmesh[i].IEN[l];
				if (l == bzmesh[i].IEN.size() - 1)
				{
					fout << "\n";
				}
				else
				{
					fout << " ";
				}
			}
			for (int j = 0; j < bzmesh[i].cmat.size(); j++)
			{
				for (int k = 0; k < bzmesh[i].cmat[j].size(); k++)
				{
					fout << bzmesh[i].cmat[j][k];
					if (k == bzmesh[i].cmat[j].size() - 1)
					{
						fout << "\n";
					}
					else
					{
						fout << " ";
					}
				}
			}

		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname1 << '\n';
	}

	// cout << "ck5" << endl;

	string fname2(fn + "bzpt.txt");
	//ofstream fout;
	fout.open(fname2.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() * 64 << "\n";
		for (int i = 0; i<bzmesh.size(); i++)
		{
			for (int j = 0; j < 64; j++)
			{
				fout << bzmesh[i].pts[j][0] << " " << bzmesh[i].pts[j][1] << " " << bzmesh[i].pts[j][2] << "\n";
			}
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname2 << '\n';
	}
}

void kernel::OutputMesh(const vector<BezierElement3D>& bzmesh, string fn, int itr)
{
	int cn[8] = { 0, 3, 15, 12, 48, 51, 63, 60 };
	string fname = fn + to_string(itr) + "_bzmesh.vtk";
	ofstream fout;
	fout.open(fname.c_str());	

	// cout << "ck2" << endl;

	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nBezier mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << 8 * bzmesh.size() << " float\n";		
		for (int i = 0; i<bzmesh.size(); i++)
		{
			for (int j = 0; j < 8; j++)
			{
				fout << bzmesh[i].pts[cn[j]][0] << " " << bzmesh[i].pts[cn[j]][1] << " " << bzmesh[i].pts[cn[j]][2] << "\n";
			}
		}
		fout << "\nCELLS " << bzmesh.size() << " " << 9 * bzmesh.size() << '\n';
		for (int i = 0; i<bzmesh.size(); i++)
		{
			fout << "8 " << 8 * i << " " << 8 * i + 1 << " " << 8 * i + 2 << " " << 8 * i + 3
				<< " " << 8 * i + 4 << " " << 8 * i + 5 << " " << 8 * i + 6 << " " << 8 * i + 7 << '\n';
		}
		fout << "\nCELL_TYPES " << bzmesh.size() << '\n';
		for (int i = 0; i<bzmesh.size(); i++)
		{
			fout << "12\n";
		}
		//fout << "POINT_DATA " << sdisp.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		//for (int i = 0; i<sdisp.size(); i++)
		//{
		//	fout << sdisp[i] << "\n";
		//}
		fout << "\nCELL_DATA " << bzmesh.size() << "\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		for (int i = 0; i < bzmesh.size(); i++)
		{
			fout << bzmesh[i].type << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
	
	// cout << "ck3" << endl;

	string fname3(fn + to_string(itr) + "_bzmeshinfo.txt");
	//ofstream fout;
	fout.open(fname3.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() << "\n";
		for (int i = 0; i<bzmesh.size(); i++)
		{
			for (int l = 0; l < bzmesh[i].IEN.size(); l++)
			{
				fout << bzmesh[i].IEN[l] + 1;
				if (l == bzmesh[i].IEN.size() - 1)
				{
					fout << "\n";
				}
				else
				{
					fout << " ";
				}
			}
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname3 << '\n';
	}

	// cout << "ck4" << endl;

	string fname1(fn + to_string(itr) + "_cmat.txt");
	//ofstream fout;
	fout.open(fname1.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() << "\n";
		for (int i = 0; i<bzmesh.size(); i++)
		{
			fout << i << " " << bzmesh[i].IEN.size() << " " << bzmesh[i].type << "\n";
			for (int l = 0; l < bzmesh[i].IEN.size(); l++)
			{
				fout << bzmesh[i].IEN[l];
				if (l == bzmesh[i].IEN.size() - 1)
				{
					fout << "\n";
				}
				else
				{
					fout << " ";
				}
			}
			for (int j = 0; j < bzmesh[i].cmat.size(); j++)
			{
				for (int k = 0; k < bzmesh[i].cmat[j].size(); k++)
				{
					fout << bzmesh[i].cmat[j][k];
					if (k == bzmesh[i].cmat[j].size() - 1)
					{
						fout << "\n";
					}
					else
					{
						fout << " ";
					}
				}
			}

		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname1 << '\n';
	}

	// cout << "ck5" << endl;

	string fname2(fn + to_string(itr) + "_bzpt.txt");
	//ofstream fout;
	fout.open(fname2.c_str());
	if (fout.is_open())
	{
		fout << bzmesh.size() * 64 << "\n";
		for (int i = 0; i<bzmesh.size(); i++)
		{
			for (int j = 0; j < 64; j++)
			{
				fout << bzmesh[i].pts[j][0] << " " << bzmesh[i].pts[j][1] << " " << bzmesh[i].pts[j][2] << "\n";
			}
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname2 << '\n';
	}
}

void kernel::run_complex_fit(string path_in)
{
	int niter(3);
	double thresh(0.25);//cube
	unsigned int i;
	double xy[3][2], nm[3], a(50.);
	TruncatedTspline_3D tt3;

	clock_t begin = clock();

	tt3.SetProblem("../ioTHS3D/hex_input/cube5");
	// tt3.SetProblem("./io/plate_input/input_CM");

	string fld("../io/kk_test/");

	string fn("base_");

	tt3.SetDomainRange(xy, nm, a);

	vector<int> dof_list(niter, 0);
	vector<double> err_list(niter, 0.);

	int itr;
	double errL2(1.e6);
	vector<BezierElement3D> bzmesh;
	for (itr = 0; itr < niter; itr++)
	{
		// vector<BezierElement3D> bzmesh;
		vector<int> IDBC;
		vector<double> gh, err;

		tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);

 		Laplace lap;
		lap.SetProblem(IDBC, gh);
		stringstream ss;
		ss << itr;
		lap.GetEqParameter(xy, nm, a);
		lap.Run(bzmesh, fld + fn + ss.str(), err);

		cout << "+++++++++++++++++++++" << endl;
		cout << err.size() << endl;
		cout << "+++++++++++++++++++++" << endl;

		errL2 = 0.;
		for (i = 0; i < err.size(); i++) errL2 += err[i];
		errL2 = sqrt(errL2);
		dof_list[itr] = IDBC.size();
		err_list[itr] = errL2;

		output_err(fld + fn + ss.str() + "_err", dof_list, err_list);

		// if (errL2 > thresh)
		if (itr > 0)
		{
			cout << itr << " refining...\n";
			//distribute error
			vector<array<double, 2>> eh(bzmesh.size());
			for (i = 0; i < bzmesh.size(); i++)
			{
				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
			}
			
			vector<array<int, 2>> rfid, gst;
			// tt3.Identify_Poisson_1(eh, err, rfid, gst);
			tt3.Identify_Laplace(eh, err, rfid, gst);			
			// tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
			// tt3.InputRefineID("../ioTHS3D/kk_test/base_" + ss.str(), rfid, gst);
			tt3.Refine(rfid, gst);
			cout << "Refining done\n";

			tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
			cout << "Output Geom done!\n";
		}
	}

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "\nElapsed time: " << elapsed_secs << "\n";

	//output error
	//output_err(fld + fn + "err", dof_list, err_list);
}

int kernel::run_neuronGrowth(string path_in, int rf_level)
{
	// int niter(4);
	int niter = 10;
	double thresh(0.25);//cube
	unsigned int i;
	double xy[3][2], nm[3], a(50.);
	TruncatedTspline_3D tt3;

	clock_t begin = clock();

	tt3.SetProblem(path_in + "controlmesh_initial");
	// tt3.SetProblem(path_in + "cube5");

	string fld(path_in);

	// string fn("base_");
	string fn("outputmesh");

	tt3.SetDomainRange(xy, nm, a);

	vector<int> dof_list(niter, 0);
	vector<double> err_list(niter, 0.);

	// cout << "ck1" << endl;
	int itr;
	double errL2(1.e6);
	vector<BezierElement3D> bzmesh;
	cout << "Start refining...\n";

	vector<int> IDBC;
	vector<double> gh;
	tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
	// tt3.GetBezierMesh(bzmesh);

	// vector<BezierElement3D> bzmesh_old = bzmesh;
	vector<BezierElement3D> bzmesh_old;

	// vector<double> phi = readVectorFromFile("../ioTHS3D/phi.txt", false);
	vector<double> phi = readVectorFromFile(path_in + "phi.txt", false);
	// int sum_of_elems = accumulate(phi.begin(), phi.end(),
        //                         decltype(phi)::value_type(0));
	// cout << "#refine phi read: " << sum_of_elems << endl;
	// vector<double> phi_old = phi;
	vector<double> phi_old;

	// for (itr = 0; itr <= niter; itr++)
	itr = 0;
	while (tt3.getLevels() <= rf_level)
	{
		if (itr >= niter) {
			cerr << "Something went wrong, refine looping too many times." << endl;
		}

		cout << "+++++++++++++++++++++" << endl;
		cout << "Refine iter " << itr << "| Current level: " << tt3.getLevels() << "...\n";

		// vector<BezierElement3D> bzmesh_old = bzmesh;
		bzmesh_old = bzmesh;
		// vector<double> phi_old = phi;
		phi_old = phi;
		// vector<BezierElement3D> bzmesh;
		// vector<int> IDBC;
		// // vector<double> gh, err;
		// vector<double> gh, err(bzmesh.size(), 0);
		vector<double> err(bzmesh.size(), 0);

		// cout << "level: " << tt3.getLevels() << endl;
		if(itr > 0) {
			cout << "Reading bzmesh...\n";
			tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
			// tt3.GetBezierMesh(bzmesh);
			// tt3.VisualizeControlMesh("../ioTHS3D/controlmesh");
			// OutputMesh(bzmesh, "../ioTHS3D/", itr);	
			// tt3.OutputControlPoints("../ioTHS3D/controlmesh", itr);
			// OutputMesh(bzmesh, "../ioTHS3D/");
			// tt3.OutputControlPoints("../ioTHS3D/");
			OutputMesh(bzmesh, path_in);
			tt3.OutputControlPoints(path_in);
			if (tt3.getLevels() == rf_level) {
			// if (itr == niter-1) {
				cout << niter << endl;
				return 0;		
			}		
		}
		// itr += 1;

		// int ini_bzmesh_size;
		// if (itr == 0) {
		// 	ini_bzmesh_size = bzmesh.size();
		// }

 		// Laplace lap;
		// lap.SetProblem(IDBC, gh);
		stringstream ss;
		ss << itr;
		// lap.GetEqParameter(xy, nm, a);
		// lap.Run(bzmesh, fld + fn + ss.str(), err);

		phi = InterpolateValues(bzmesh_old, phi_old, bzmesh);
		writeVectorToFile(phi, "./phi_refine.txt", false);
		// err.clear();
		err = phi;
		// cout << "ck2 " << bzmesh.size() << " " << err.size() << " " << ids.size() << endl;

		cout << "bzmesh size: " << bzmesh.size() << " phi size: " << phi.size() << endl;
		
		// lap.VisualizeError(bzmesh, err, fld + fn + ss.str());

		// cout << "+++++++++++++++++++++" << endl;
		// cout << err.size() << endl;
		// cout << "Refining iter " << itr << "...\n";
		// cout << "+++++++++++++++++++++" << endl;

		errL2 = 0.;
		for (i = 0; i < err.size(); i++) errL2 += err[i];
		errL2 = sqrt(errL2);
		dof_list[itr] = IDBC.size();
		err_list[itr] = errL2;

		cout << "err size: " << err.size() << endl;;

		// output_err(fld + fn + ss.str() + "_err", dof_list, err);

	
		// cout << itr << " refining...\n";
		// distribute error
		vector<array<double, 2>> eh(bzmesh.size());
		for (i = 0; i < bzmesh.size(); i++)
		{
			eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
		}
		
		vector<array<int, 2>> rfid, gst;
		// tt3.Identify_Poisson_1(eh, err, rfid, gst);
		tt3.Identify_Laplace(eh, err, rfid, gst);			
		// tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
		// tt3.InputRefineID(fld + fn + ss.str(), rfid, gst);
		tt3.Refine(rfid, gst);
		// cout << "Refining done\n";

		tt3.OutputGeom_All(fld + fn + "_" + ss.str() + "_geom");
		// cout << "Output Geom done!\n";

		cout << "Refine iter " << itr << " done!\n";

		// cout << "ck0" << endl;

		// if (itr == niter) { // update bzmesh for outputmesh
		// 	tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
		// }
		// tt3.VisualizeControlMesh("../ioTHS3D/controlmesh");
		// tt3.OutputCM(itr, "../ioTHS3D/controlmesh");

		// OutputMesh(bzmesh, "../ioTHS3D/");
		// cout << eh.size() << " " << err.size() << endl;

		// tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
		// // tt3.OutputCM(itr, "../ioTHS3D/controlmesh");
		// OutputMesh(bzmesh, "../ioTHS3D/", itr);	
		// tt3.OutputControlPoints("../ioTHS3D/controlmesh", itr);

		itr++;
	}
	
	// tt3.OutputControlPoints("../ioTHS3D/controlmesh");
	// tt3.VisualizeControlMesh("../ioTHS3D/controlmesh");

	// tt3.OutputCM_allLevel("../ioTHS3D/controlmesh");
	// cout << "Writing Hierarachical mesh ..." << endl;
	// tt3.VisualizeControlMesh_hierarchical("../ioTHS3D/controlmesh");
	// cout << "Writing Tmesh ..." << endl;
	// tt3.VisualizeTMesh("../ioTHS3D/controlmesh");
	// cout << "Writing CM ..." << endl;
	// tt3.OutputCM("../ioTHS3D/controlmesh");

	// cout << "ck1" << endl;
	// tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
	// OutputMesh(bzmesh, "../ioTHS3D/");

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "\nElapsed time: " << elapsed_secs << "\n";

	return 0;
}

// A small struct for convenience
struct NeighborInfo {
    int index;      // index in the old mesh
    double distance;
};

// This function returns a sorted list (by ascending distance) of the n nearest neighbors
vector<NeighborInfo> FindKNearestNeighbors(const vector<BezierElement3D>& bzmesh_old,
												const BezierElement3D& target_element,
												int K)
{
    vector<NeighborInfo> neighbors;
    neighbors.reserve(bzmesh_old.size());

    // Compute distance from target_element to each old element
    for (int i = 0; i < (int)bzmesh_old.size(); ++i) {
        double distance = 0.0;
        for (int j = 0; j < 3; ++j) {
            double diff = target_element.pts[0][j] - bzmesh_old[i].pts[0][j];
            distance += diff * diff;
        }
        // distance = sqrt(distance);

        NeighborInfo info;
        info.index = i;
        info.distance = distance;
        neighbors.push_back(info);
    }

    // Sort by distance ascending
    sort(neighbors.begin(), neighbors.end(),
              [](const NeighborInfo& a, const NeighborInfo& b) {
                  return a.distance < b.distance;
              });

    // Truncate to the first K neighbors if the array is larger
    if ((int)neighbors.size() > K) {
        neighbors.resize(K);
    }

    return neighbors;
}

// // Function to find the index of the nearest neighbor in the old mesh
// int kernel::FindNearestNeighbor(const vector<BezierElement3D>& bzmesh_old, const BezierElement3D& target_element, double& min_distance, const vector<double>& phi_old)
// {
// 	int nearest_index = 0;
// 	min_distance = numeric_limits<double>::max();
// 	// double min_distance = numeric_limits<double>::max();

// 	for (int i = 0; i < bzmesh_old.size(); ++i) {
// 		// if (phi_old[i] != 0) {
// 			// Calculate the Euclidean distance between target_element and elements in bzmesh_old
// 			double distance = 0.0;
// 			for (int j = 0; j < 3; ++j)
// 			{
// 				double diff = target_element.pts[0][j] - bzmesh_old[i].pts[0][j];
// 				distance += diff * diff;
// 			}
// 			distance = sqrt(distance);

// 			// Update nearest neighbor if a closer one is found
// 			if (distance < min_distance) {
// 				min_distance = distance;
// 				nearest_index = i;
// 			}
// 		// }
// 	}

// 	return nearest_index;
// }

vector<double> kernel::InterpolateValues(const vector<BezierElement3D>& bzmesh_old,
                                         const vector<double>& phi_old,
                                         const vector<BezierElement3D>& bzmesh_new)
{
	int K = 6;
	double DIST_THRESHOLD = 12.0;

    // Ensure input sizes match
    assert(bzmesh_old.size() == phi_old.size());

    // Initialize output
    vector<double> phi_new(bzmesh_new.size(), 0.0);

    // For each new element, find the K nearest neighbors
    for (int i = 0; i < (int)bzmesh_new.size(); ++i) {
        auto neighbors = FindKNearestNeighbors(bzmesh_old, bzmesh_new[i], K);

        // Check if any of the neighbors has a non-zero value
        for (const auto& nb : neighbors) {
            if (nb.distance <= DIST_THRESHOLD && phi_old[nb.index] != 0.0) {
                phi_new[i] = 1.0; // Set phi_new[i] to 1 if any neighbor is non-zero
                break;           // No need to check further neighbors
            }
        }
    }

    return phi_new;
}

// vector<double> kernel::InterpolateValues(const vector<BezierElement3D>& bzmesh_old,
//                                          const vector<double>& phi_old,
//                                          const vector<BezierElement3D>& bzmesh_new)
// {
//     // 1) Choose how many neighbors (K) you want to average
//     const int K = 6;               // e.g. 6 nearest neighbors
//     const double DIST_THRESHOLD = 12.0;

//     // 2) Initialize output
//     vector<double> phi_new(bzmesh_new.size(), 0.0);
//     vector<bool>   updated(bzmesh_new.size(), false);

//     // 3) For each new element, find the K nearest neighbors from the old mesh
//     for (int i = 0; i < (int)bzmesh_new.size(); ++i) {
//         // Gather the K nearest neighbors
//         auto neighbors = FindKNearestNeighbors(bzmesh_old, bzmesh_new[i], K);

//         // If the closest neighbor is beyond the threshold, skip
//         if (!neighbors.empty()) {
//             double closestDist = neighbors[0].distance;
//             if (!updated[i] && closestDist <= DIST_THRESHOLD) {
//                 // 4) Compute the weighted average phi of these neighbors
//                 double weighted_sum_phi = 0.0;
//                 double weight_sum = 0.0;

//                 for (auto& nb : neighbors) {
//                     // Only consider neighbors within the distance threshold
//                     if (nb.distance <= DIST_THRESHOLD) {
//                         double weight = 1.0 / (nb.distance + 1e-6); // Avoid division by zero
//                         weighted_sum_phi += phi_old[nb.index] * weight;
//                         weight_sum += weight;
//                     }
//                 }

//                 if (weight_sum > 0) {
//                     phi_new[i] = weighted_sum_phi / weight_sum; // Weighted average
//                     updated[i] = true;
//                 }
//             }
//         }
//     }

//     return phi_new;
// }

// vector<double> kernel::InterpolateValues(const vector<BezierElement3D>& bzmesh_old,
//                                      const vector<double>& phi_old,
//                                      const vector<BezierElement3D>& bzmesh_new)
// {
//     // Initialize phi_new with 0.0
//     vector<double> phi_new(bzmesh_new.size(), 0.0);

//     // Track updated indices
//     vector<bool> updated(bzmesh_new.size(), false);

//     for (int i = 0; i < bzmesh_new.size(); ++i) {
//         // Find the nearest neighbor in the old mesh for each element in bzmesh_new
//         double dist(10);
//         int nearest_index = FindNearestNeighbor(bzmesh_old, bzmesh_new[i], dist, phi_old);

//         // Update phi_new only if it has not been updated and distance is within the threshold
//         if (!updated[i] && dist <= 6) {
//             phi_new[i] = phi_old[nearest_index];
//             updated[i] = true; // Mark this index as updated
//         }
//     }

//     // No need to explicitly set unupdated elements to 0 as phi_new is already initialized to 0.0

//     return phi_new;
// }

// // Function to perform interpolation from old mesh to new mesh
// vector<double> kernel::InterpolateValues(const vector<BezierElement3D>& bzmesh_old,
//                                      const vector<double>& phi_old,
//                                      const vector<BezierElement3D>& bzmesh_new)
// {
// 	vector<double> phi_new(bzmesh_new.size(), 0.0);

// 	for (int i = 0; i < bzmesh_new.size(); ++i) {
// 		// Find the nearest neighbor in the old mesh for each element in bzmesh_new
// 		double dist(10);
// 		int nearest_index = FindNearestNeighbor(bzmesh_old, bzmesh_new[i], dist, phi_old);
// 		// cout << dist << endl;
// 		// Interpolate the value based on the nearest neighbor
// 		// phi_new[i] = phi_old[nearest_index];
// 		if (dist <= 6) {
// 			phi_new[i] = phi_old[nearest_index];
// 		} else {
// 			phi_new[i] = 0;
// 		}
			
// 	}

// 	return phi_new;
// }

void kernel::writeVectorToFile(const vector<double>& data, const string& filename, bool binary) {
	ofstream outfile;

	if (binary) {
		outfile.open(filename, ios::out | ios::binary);
	} else {
		outfile.open(filename);
	}

	if (!outfile) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}

	if (binary) {
		outfile.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(double));
	} else {
		for (const auto& value : data) {
			outfile << value << " ";
		}
	}

	cout << "Vector successfully written to " << filename << endl;
	outfile.close();
}

vector<double> kernel::readVectorFromFile(const string& filename, bool binary) {
	ifstream infile;

	if (binary) {
		infile.open(filename, ios::in | ios::binary);
	} else {
		infile.open(filename);
	}

	if (!infile) {
		cerr << "Error opening file: " << filename << endl;
		return {};
	}

	vector<double> data;

	if (binary) {
		infile.seekg(0, ios::end);
		size_t fileSize = infile.tellg();
		infile.seekg(0, ios::beg);

		data.resize(fileSize / sizeof(double));
		infile.read(reinterpret_cast<char*>(data.data()), fileSize);
	} else {
		double value;

		while (infile >> value) {
			data.push_back(value);
		}
	}

	cout << "Vector successfully read from " << filename << endl;
	infile.close();

	return data;
}

// void kernel::run_complex_glb()
// {
// 	int niter(1);
// 	unsigned int i;
// 	double xy[3][2], nm[3], a(50.);
// 	TruncatedTspline_3D tt3;
// 	//tt3.InitializeMesh("../io/hex_input/cube_coarse_0");
// 	tt3.InitializeMesh("../io/hex_input/rod");
// 	//tt3.InitializeMesh("../io/hex_input/hook");
// 	//tt3.InitializeMesh("../io/hex_input/cad1");
// 	//tt3.InitializeMesh("../io/hex_input/statue");
// 	//tt3.InitializeMesh("../io/hex_input/gear");

// 	//tt3.SetProblem("../io/hex_input/cube_dense");
// 	//tt3.SetProblem("../io/hex_input/cube3");
// 	//tt3.SetProblem("../io/hex_input/gear");

// 	string fld("../io/global3/");
// 	//string fld("../io/benchmark3/");

// 	//string fn("cube0_4");
// 	string fn("rod1_2");
// 	//string fn("hook1_");
// 	//string fn("base0_");
// 	//string fn("head1_");
// 	//string fn("gear1_");

// 	tt3.run_Global(2, fld + fn);
// 	cout << "Done global refinement!\n";
// 	//getchar();

// 	//double remove[3][2];
// 	//tt3.GetRemoveRegion(remove);
// 	tt3.SetDomainRange(xy, nm, a);

// 	//tt3.OutputRemoveCM(fld+fn,remove);
// 	//cout << "done\n";
// 	//getchar();

// 	vector<int> dof_list(niter, 0);
// 	vector<double> h_max(niter, 0.);
// 	vector<double> err_list(niter, 0.);
// 	for (int itr = 0; itr < niter; itr++)
// 	{
// 		vector<BezierElement3D> bzmesh;
// 		vector<int> IDBC;
// 		vector<double> gh, err;
// 		tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);

// 		Laplace lap;
// 		lap.SetProblem(IDBC, gh);
// 		stringstream ss;
// 		ss << itr;
// 		//lap.GetRemoveRegion(remove);
// 		lap.GetEqParameter(xy, nm, a);
// 		lap.Run(bzmesh, fld + fn, err);

// 		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());
// 		//cout << "Done visualize Bezier!\n";
// 		//getchar();

// 		double errL2(0.);
// 		for (i = 0; i < err.size(); i++) errL2 += err[i];
// 		errL2 = sqrt(errL2);
// 		dof_list[itr] = IDBC.size();
// 		h_max[itr] = tt3.MaxElementSize(bzmesh);
// 		err_list[itr] = errL2;
// 		//cout << "DOF: " << IDBC.size() << "\n";
// 		//cout << "L2-norm error: " << errL2 << "\n";
// 		//getchar();

// 		output_err(fld + fn + "_err", dof_list, err_list);
// 		output_err(fld + fn + "_err_h", h_max, err_list);
// 	}
// }

// void kernel::run_Bspline()
// {
// 	int niter(1);
// 	unsigned int i;
// 	string fld("../io/BSP/");
// 	//string fn("ZBSP3_");
// 	string fn("ZBSP3_");
// 	TruncatedTspline_3D tt3;
// 	vector<BezierElement3D> bzmesh;
// 	vector<int> IDBC;
// 	vector<double> gh, err;
// 	int nsp(4);
// 	int nex[3] = { nsp, nsp, nsp };
// 	//tt3.CreateBsplines(bzmesh, IDBC, gh);
// 	tt3.CreateBsplines(nex);
// 	int nref(3);
// 	for (int i = 0; i < nref; i++)
// 	{
// 		tt3.Bsplines_Refine();
// 	}
// 	tt3.Bspline_BezierExtract(bzmesh, IDBC, gh);

// 	//tt3.VisualizeBezier(bzmesh, fld + fn);
// 	//cout << "Done visualize Bezier!\n";
// 	//getchar();

// 	Laplace lap;
// 	lap.SetProblem(IDBC, gh);
// 	lap.Run(bzmesh, fld + fn, err);

// 	vector<int> dof_list(niter, 0);
// 	vector<double> h_max(niter, 0.);
// 	vector<double> err_list(niter, 0.);
// 	double errL2(0.);
// 	for (i = 0; i < err.size(); i++) errL2 += err[i];
// 	errL2 = sqrt(errL2);
// 	dof_list[0] = IDBC.size();
// 	h_max[0] = tt3.MaxElementSize(bzmesh);
// 	err_list[0] = errL2;

// 	output_err(fld + fn + "_err", dof_list, err_list);
// 	output_err(fld + fn + "_err_h", h_max, err_list);
// }

// void kernel::run_Bspline_fit()
// {
// 	int niter(1);
// 	unsigned int i;
// 	double xy[3][2], nm[3], a(50.);
// 	string fld("../io/BSP/");
// 	string fn("BSP2_fit4");
// 	TruncatedTspline_3D tt3;
// 	vector<BezierElement3D> bzmesh;
// 	vector<int> IDBC;
// 	vector<double> gh, err;
// 	int nsp(2);
// 	int nex[3] = { nsp, nsp, nsp };
// 	//tt3.CreateBsplines(bzmesh, IDBC, gh);
// 	tt3.CreateBsplines(nex);
// 	int nref(4);
// 	for (int i = 0; i < nref; i++)
// 	{
// 		tt3.Bsplines_Refine();
// 	}
// 	tt3.SetDomainRange_Bsplines(xy, nm, a);//for solution 7
// 	tt3.Bspline_BezierExtract_fit(bzmesh, IDBC, gh);
// 	tt3.FittingBC_Bsplines(bzmesh, IDBC, gh);

// 	//tt3.VisualizeBezier(bzmesh, fld + fn);
// 	//cout << "Done visualize Bezier!\n";
// 	//getchar();

// 	Laplace lap;
// 	lap.SetProblem(IDBC, gh);
// 	lap.GetEqParameter(xy, nm, a);//for solution 7
// 	lap.Run(bzmesh, fld + fn, err);

// 	vector<int> dof_list(niter, 0);
// 	vector<double> h_max(niter, 0.);
// 	vector<double> err_list(niter, 0.);
// 	double errL2(0.);
// 	for (i = 0; i < err.size(); i++) errL2 += err[i];
// 	errL2 = sqrt(errL2);
// 	dof_list[0] = IDBC.size();
// 	h_max[0] = tt3.MaxElementSize(bzmesh);
// 	err_list[0] = errL2;

// 	output_err(fld + fn + "_err", dof_list, err_list);
// 	output_err(fld + fn + "_err_h", h_max, err_list);
// }

// void kernel::run_benchmark()
// {
// 	int niter(4);
// 	unsigned int i;
// 	TruncatedTspline_3D tt3;
// 	//tt3.SetProblem("../io/hex_input/cube_dense");
// 	//tt3.SetProblem("../io/hex_input/cube9");
// 	//tt3.SetProblem("../io/hex_input/cube_coarse_0");
// 	tt3.SetProblem("../io/hex_input/cube_coarse_2");
// 	string fld("../io/benchmark4/");
// 	string fn("cube_THS1_");
// 	//string fn("cube_HS1_");

// 	//double remove[3][2];
// 	//tt3.GetRemoveRegion(remove);
// 	//tt3.OutputRemoveCM(fld+fn,remove);
// 	//cout << "done\n";
// 	//getchar();

// 	vector<int> dof_list(niter, 0);
// 	vector<double> err_list(niter, 0.);
// 	vector<array<int, 2>> ebc, pbc;
// 	vector<double> edisp, pdisp;
// 	//tt3.SetInitialBC(ebc, edisp, pbc, pdisp);
// 	for (int itr = 0; itr < niter; itr++)
// 	{
// 		vector<BezierElement3D> bzmesh;
// 		vector<int> IDBC;
// 		vector<double> gh, err;
// 		//tt3.SetBC(ebc, edisp, pbc, pdisp);
// 		//cout << pbc.size() << "\n";
// 		//cout << "interface...\n";
// 		tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
// 		//tt3.AnalysisInterface_Laplace(pbc,pdisp,bzmesh, IDBC, gh);
// 		//cout << "interface done!\n";
// 		//getchar();

// 		//cout << "npt: " << IDBC.size() << "\n";
// 		//getchar();

// 		Laplace lap;
// 		lap.SetProblem(IDBC, gh);
// 		stringstream ss;
// 		ss << itr;
// 		//cout << "before simulation...\n";
// 		//getchar();
// 		//lap.GetRemoveRegion(remove);
// 		lap.Run(bzmesh, fld + fn + ss.str(), err);
// 		//cout << "simulation done!\n";
// 		//getchar();

// 		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

// 		double errL2(0.);
// 		for (i = 0; i < err.size(); i++) errL2 += err[i];
// 		errL2 = sqrt(errL2);
// 		dof_list[itr] = IDBC.size();
// 		err_list[itr] = errL2;
// 		//cout << "DOF: " << IDBC.size() << "\n";
// 		//cout << "L2-norm error: " << errL2 << "\n";
// 		//getchar();

// 		output_err(fld + fn + ss.str() + "_err", dof_list, err_list);

// 		if (itr < niter - 1)
// 		{
// 			cout << itr << " refining...\n";
// 			//distribute error
// 			vector<array<double, 2>> eh(bzmesh.size());
// 			for (i = 0; i < bzmesh.size(); i++)
// 			{
// 				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
// 			}
// 			vector<array<int, 2>> rfid, gst;
// 			//tt3.Identify_Poisson_1(eh, err, rfid, gst);
// 			//tt3.Identify_Laplace(eh, err, rfid, gst);
// 			//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
// 			tt3.InputRefineID("../io/benchmark3/cube_THS1_" + ss.str(), rfid, gst);
// 			tt3.Refine(rfid, gst);
// 			cout << "Refining done\n";
// 			//tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
// 			//cout << "Output Geom done!\n";
// 			//getchar();
// 		}
// 	}

// 	//output error
// 	//output_err(fld + fn + "err", dof_list, err_list);
// }

// void kernel::run_leastsquare()
// {
// 	int niter(2);
// 	unsigned int i;
// 	TruncatedTspline_3D tt3;
// 	//tt3.SetProblem("../io/hex_input/cube_dense");
// 	//tt3.SetProblem("../io/hex_input/cube4");
// 	tt3.SetProblem("../io/hex_input/cube_coarse_2");
// 	string fld("../io/least_square1/");
// 	//string fn("cube_THS2_2");
// 	string fn("cube_HS2_2");

// 	//tt3.GlobalRefine(2);

// 	for (int itr = 0; itr < niter; itr++)
// 	{
// 		cout << itr << " refining...\n";
// 		vector<array<int, 2>> rfid, gst;
// 		//tt3.Identify_Poisson_1(eh, err, rfid, gst);
// 		//tt3.Identify_Laplace(eh, err, rfid, gst);
// 		//tt3.Identify_LeastSquare(rfid, gst);
// 		tt3.Identify_LeastSquare_Line(rfid, gst);
// 		//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
// 		//tt3.InputRefineID("../io/benchmark3/cube_THS1_" + ss.str(), rfid, gst);
// 		tt3.Refine(rfid, gst);
// 		stringstream ss;
// 		ss << itr;
// 		cout << "Refining done\n";
// 		//tt3.OutputGeom_All(fld + fn + ss.str());
// 		//cout << "Output Geom done!\n";
// 		//getchar();
// 	}

// 	//cout << "Refining done\n";
// 	//tt3.OutputGeom_All(fld + fn + "_geom");
// 	//cout << "geom done!\n";
// 	//getchar();

// 	vector<BezierElement3D> bzmesh;
// 	vector<int> IDBC;
// 	vector<double> gh, err;
// 	//tt3.AnalysisInterface_LeastSquare(bzmesh, IDBC, gh);
// 	tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
// 	LeastSquare lap;
// 	lap.SetProblem(IDBC, gh);
// 	lap.Run(bzmesh, fld + fn, err);

// 	//output error
// 	//output_err(fld + fn + "err", dof_list, err_list);
// }

// void kernel::run_leastsquare_1()
// {
// 	int niter(4);
// 	unsigned int i;
// 	TruncatedTspline_3D tt3;
// 	//tt3.SetProblem("../io/hex_input/cube_dense");
// 	//tt3.SetProblem("../io/hex_input/cube9");
// 	//tt3.SetProblem("../io/hex_input/cube_coarse_0");
// 	tt3.SetProblem("../io/hex_input/cube_coarse_2");
// 	string fld("../io/least_square1/");
// 	string fn("cube_THS1_");
// 	//string fn("cube_HS1_");

// 	//double remove[3][2];
// 	//tt3.GetRemoveRegion(remove);
// 	//tt3.OutputRemoveCM(fld+fn,remove);
// 	//cout << "done\n";
// 	//getchar();

// 	vector<int> dof_list(niter, 0);
// 	vector<double> err_list(niter, 0.);
// 	vector<array<int, 2>> ebc, pbc;
// 	vector<double> edisp, pdisp;
// 	//tt3.SetInitialBC(ebc, edisp, pbc, pdisp);
// 	for (int itr = 0; itr < niter; itr++)
// 	{
// 		vector<BezierElement3D> bzmesh;
// 		vector<int> IDBC;
// 		vector<double> gh, err;
// 		//tt3.SetBC(ebc, edisp, pbc, pdisp);
// 		//cout << pbc.size() << "\n";
// 		//cout << "interface...\n";
// 		//tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
// 		//tt3.AnalysisInterface_Laplace(pbc,pdisp,bzmesh, IDBC, gh);
// 		tt3.AnalysisInterface_LeastSquare(bzmesh, IDBC, gh);
// 		//cout << "interface done!\n";
// 		//getchar();

// 		//cout << "npt: " << IDBC.size() << "\n";
// 		//getchar();

// 		Laplace lap;
// 		lap.SetProblem(IDBC, gh);
// 		stringstream ss;
// 		ss << itr;
// 		//cout << "before simulation...\n";
// 		//getchar();
// 		//lap.GetRemoveRegion(remove);
// 		lap.Run(bzmesh, fld + fn + ss.str(), err);
// 		//cout << "simulation done!\n";
// 		//getchar();

// 		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

// 		double errL2(0.);
// 		for (i = 0; i < err.size(); i++) errL2 += err[i];
// 		errL2 = sqrt(errL2);
// 		dof_list[itr] = IDBC.size();
// 		err_list[itr] = errL2;
// 		//cout << "DOF: " << IDBC.size() << "\n";
// 		//cout << "L2-norm error: " << errL2 << "\n";
// 		//getchar();

// 		output_err(fld + fn + ss.str() + "_err", dof_list, err_list);

// 		if (itr < niter - 1)
// 		{
// 			cout << itr << " refining...\n";
// 			//distribute error
// 			vector<array<double, 2>> eh(bzmesh.size());
// 			for (i = 0; i < bzmesh.size(); i++)
// 			{
// 				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
// 			}
// 			vector<array<int, 2>> rfid, gst;
// 			//tt3.Identify_Poisson_1(eh, err, rfid, gst);
// 			//tt3.Identify_Laplace(eh, err, rfid, gst);
// 			//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
// 			tt3.InputRefineID("../io/benchmark3/cube_THS1_" + ss.str(), rfid, gst);
// 			tt3.Refine(rfid, gst);
// 			cout << "Refining done\n";
// 			//tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
// 			//cout << "Output Geom done!\n";
// 			//getchar();
// 		}
// 	}

// 	//output error
// 	//output_err(fld + fn + "err", dof_list, err_list);
// }

// void kernel::run_pipeline()
// {
// 	int niter(5);
// 	unsigned int i;
// 	TruncatedTspline_3D tt3;
// 	tt3.SetProblem("../io/pipeline/neuron");
// 	string fld("../io/pipeline/");
// 	string fn("neuron1");

// 	//double remove[3][2];
// 	//tt3.GetRemoveRegion(remove);
// 	//tt3.OutputRemoveCM(fld+fn,remove);
// 	//cout << "done\n";
// 	//getchar();

// 	vector<BezierElement3D> bzmesh;
// 	tt3.PipelineBezierExtract(bzmesh);

// 	Laplace lap;
// 	lap.PipelineTmp(bzmesh,fld+fn);
// }

void kernel::output_err(string fn, const vector<int>& dof, const vector<double>& err)
{
	string fname = fn + ".txt";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		// for (i = 0; i < dof.size(); i++)
		// {
		// 	fout << dof[i] << " " << err[i] << "\n";
		// }
		for (i = 0; i < err.size(); i++)
		{
			fout << err[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void kernel::output_err(string fn, const vector<double>& dof, const vector<double>& err)
{
	string fname = fn + ".txt";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		for (i = 0; i < dof.size(); i++)
		{
			fout << dof[i] << " " << err[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}