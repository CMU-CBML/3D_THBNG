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

	// std::cout << "ck2" << std::endl;

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
	
	// std::cout << "ck3" << std::endl;

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

	// std::cout << "ck4" << std::endl;

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

	// std::cout << "ck5" << std::endl;

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

	// std::cout << "ck2" << std::endl;

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
	
	// std::cout << "ck3" << std::endl;

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

	// std::cout << "ck4" << std::endl;

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

	// std::cout << "ck5" << std::endl;

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

		std::cout << "+++++++++++++++++++++" << std::endl;
		std::cout << err.size() << std::endl;
		std::cout << "+++++++++++++++++++++" << std::endl;

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


int kernel::run_neuronGrowth_optionA(const string& path_in, int level)
{
    // 1) Instantiate T-spline object and load initial mesh
    TruncatedTspline_3D tt3;
    tt3.SetProblem(path_in + "controlmesh_initial");

    // 2) Read the initial 'phi.txt' which indicates elements to refine
    vector<double> phi_old = readVectorFromFile(path_in + "phi.txt", false);

    // 3) Prepare containers for the current mesh
    vector<BezierElement3D> bzmesh;
    vector<int> IDBC;         // For boundary conditions
    vector<double> gh;        // For geometry/function values if needed

    // 4) Build the initial mesh representation
    tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);

    // Keep a copy of the old mesh before refinement
    vector<BezierElement3D> bzmesh_old = bzmesh; 

    // We want 3 local refinements total
    const int numRefinements = level; 

    for (int pass = 0; pass < numRefinements; pass++) {
        cout << "\n===== Local Refinement Pass " << pass+1 
             << " of " << numRefinements << " =====\n";

        // (A) Interpolate old phi onto the current mesh
        //     So that we know which elements in the new mesh should be flagged
        vector<double> phi_new = InterpolateValues(bzmesh_old, phi_old, bzmesh);

        // (B) Identify which elements to refine using 'Identify_Laplace'
        //     We'll build 'eh' to map each element in 'bzmesh' to a level/index
        vector<array<double, 2>> eh(bzmesh.size());
        for (size_t i = 0; i < bzmesh.size(); i++) {
            // 'prt' typically stores hierarchical or indexing info in [0], [1]
            eh[i][0] = bzmesh[i].prt[0];
            eh[i][1] = bzmesh[i].prt[1];
        }

        // We'll use 'phi_new' as the error array
        vector<array<int, 2>> rfid, gst;
        tt3.Identify_Laplace(eh, phi_new, rfid, gst);

        // (C) Refine the mesh based on 'rfid'/'gst' 
        tt3.Refine(rfid, gst);

        // (D) Re-generate the updated mesh with 'AnalysisInterface_Poisson_1'
        //     This step updates 'bzmesh' so we see the newly refined elements
        bzmesh_old = bzmesh;      // Keep a copy for the next pass interpolation
        tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);

        // (E) Now we want to preserve the newly computed phi for the next pass
        //     So that we can re-check. We'll store 'phi_new' in 'phi_old'
        phi_old = phi_new;
    }

    // (F) (Optional) Output the final mesh for visualization
    tt3.OutputGeom_All(path_in + "final_geom");

    cout << "All " << numRefinements << " local refinement passes complete.\n";

    return 0;
}

int kernel::run_neuronGrowth(string path_in, int level)
{
	// int niter(2);
	int niter = level;
	// double thresh(0.25);//cube
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

	// std::cout << "ck1" << std::endl;
	int itr;
	double errL2(1.e6);
	vector<BezierElement3D> bzmesh;
	cout << "Start refining...\n";

	vector<int> IDBC;
	vector<double> gh;
	tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
	// tt3.GetBezierMesh(bzmesh);

	vector<BezierElement3D> bzmesh_old = bzmesh;

	// vector<double> phi = readVectorFromFile("../ioTHS3D/phi.txt", false);
	vector<double> phi = readVectorFromFile(path_in + "phi.txt", false);
	// int sum_of_elems = std::accumulate(phi.begin(), phi.end(),
        //                         decltype(phi)::value_type(0));
	// std::cout << "#refine phi read: " << sum_of_elems << std::endl;
	vector<double> phi_old = phi;

	for (itr = 0; itr <= niter; itr++)
	// itr = 0;
	// while (tt3.getLevels() < 2)
	{
		std::cout << "+++++++++++++++++++++" << std::endl;
		cout << "Refine iter " << itr << "...\n";

		// vector<BezierElement3D> bzmesh_old = bzmesh;
		// vector<double> phi_old = phi;
		// vector<BezierElement3D> bzmesh;
		// vector<int> IDBC;
		// // vector<double> gh, err;
		// vector<double> gh, err(bzmesh.size(), 0);
		vector<double> err(bzmesh.size(), 0);

		// std::cout << "level: " << tt3.getLevels() << std::endl;
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
			// if (tt3.getLevels() == 3) {
			if (itr == niter-1) {
				std::cout << niter << std::endl;
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
		writeVectorToFile(phi, path_in + "phi_refine.txt", false);
		// err.clear();
		err = phi;
		// std::cout << "ck2 " << bzmesh.size() << " " << err.size() << " " << ids.size() << std::endl;

		std::cout << "bzmesh size: " << bzmesh.size() << " phi size: " << phi.size() << std::endl;
		
		// lap.VisualizeError(bzmesh, err, fld + fn + ss.str());

		// std::cout << "+++++++++++++++++++++" << std::endl;
		// std::cout << err.size() << std::endl;
		// cout << "Refining iter " << itr << "...\n";
		// std::cout << "+++++++++++++++++++++" << std::endl;

		errL2 = 0.;
		for (i = 0; i < err.size(); i++) errL2 += err[i];
		errL2 = sqrt(errL2);
		dof_list[itr] = IDBC.size();
		err_list[itr] = errL2;

		std::cout << "err size: " << err.size() << std::endl;;

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

		// std::cout << "ck0" << std::endl;

		// if (itr == niter) { // update bzmesh for outputmesh
		// 	tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
		// }
		// tt3.VisualizeControlMesh("../ioTHS3D/controlmesh");
		// tt3.OutputCM(itr, "../ioTHS3D/controlmesh");

		// OutputMesh(bzmesh, "../ioTHS3D/");
		// std::cout << eh.size() << " " << err.size() << std::endl;

		// tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
		// // tt3.OutputCM(itr, "../ioTHS3D/controlmesh");
		// OutputMesh(bzmesh, "../ioTHS3D/", itr);	
		// tt3.OutputControlPoints("../ioTHS3D/controlmesh", itr);

	}
	// tt3.OutputControlPoints("../ioTHS3D/controlmesh");
	// tt3.VisualizeControlMesh("../ioTHS3D/controlmesh");

	// tt3.OutputCM_allLevel("../ioTHS3D/controlmesh");
	// std::cout << "Writing Hierarachical mesh ..." << std::endl;
	// tt3.VisualizeControlMesh_hierarchical("../ioTHS3D/controlmesh");
	// std::cout << "Writing Tmesh ..." << std::endl;
	// tt3.VisualizeTMesh("../ioTHS3D/controlmesh");
	// std::cout << "Writing CM ..." << std::endl;
	// tt3.OutputCM("../ioTHS3D/controlmesh");

	// std::cout << "ck1" << std::endl;
	// tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
	// OutputMesh(bzmesh, "../ioTHS3D/");

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "\nElapsed time: " << elapsed_secs << "\n";

	return 0;
}

// int kernel::run_neuronGrowth(string path_in, int level)
// {
//     // 1) Create and initialize the T-spline object
//     TruncatedTspline_3D tt3;
//     tt3.SetProblem(path_in + "controlmesh_initial");

//     // 2) Load or compute the initial Bezier mesh
//     std::vector<BezierElement3D> bzmesh;
//     std::vector<int> IDBC;
//     std::vector<double> gh;
//     tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);

//     // 3) Read 'phi.txt' indicating which elements are initially flagged
//     //    for local refinement
//     std::vector<double> phi_old = readVectorFromFile(path_in + "phi.txt", false);
//     std::cout << "Vector successfully read from " << path_in << "phi.txt"
//               << " with size: " << phi_old.size() << "\n";

//     // We'll iterate 'level' times for local refinement
//     // For each pass, we refine the mesh and re-check the new mesh structure
//     for (int pass = 0; pass < level; pass++) {
//         std::cout << "\n===== Local Refinement Pass " 
//                   << (pass + 1) << " of " << level << " =====\n";

//         // (A) Interpolate old phi onto the current 'bzmesh'
//         //     (Old mesh is effectively 'bzmesh' from previous iteration, so
//         //      we keep 'bzmesh' consistent or store a copy if needed).
//         //     If this is the first pass, 'bzmesh' is the original mesh.
//         //     Otherwise, it’s the newly refined mesh from the last pass.
//         std::vector<double> phi_new = InterpolateValues(
//             bzmesh,  // old mesh (or the same if first pass)
//             phi_old, // old phi
//             bzmesh   // new mesh is also 'bzmesh' if we re-check
//         );

//         // (B) Build 'eh' to map each element in 'bzmesh' to hierarchical indexes
//         std::vector<std::array<double, 2>> eh(bzmesh.size());
//         for (size_t i = 0; i < bzmesh.size(); i++) {
//             eh[i][0] = bzmesh[i].prt[0];
//             eh[i][1] = bzmesh[i].prt[1];
//         }

//         // (C) Identify which elements to refine using 'phi_new' as the error metric
//         std::vector<std::array<int, 2>> rfid, gst;
//         tt3.Identify_Laplace(eh, phi_new, rfid, gst);

//         // (D) Refine the flagged elements
//         if (!rfid.empty() || !gst.empty()) {
//             tt3.Refine(rfid, gst);
//             std::cout << "Refining done. rfid.size()=" << rfid.size() 
//                       << ", gst.size()=" << gst.size() << "\n";
//         } else {
//             std::cout << "No elements flagged for refinement this pass.\n";
//         }

//         // (E) Output geometry or mesh if desired (for debugging or final output)
//         //     'pass' indicates which iteration. 
//         //     E.g., "geom_pass0", "geom_pass1", "geom_pass2", etc.
//         std::stringstream ss;
//         ss << "geom_pass" << pass;
//         tt3.OutputGeom_All(path_in + ss.str());

//         // (F) Re-generate the updated Bezier mesh to reflect the refined structure
//         //     so next pass re-checks a newly refined mesh.
//         bzmesh.clear(); // optional, since 'AnalysisInterface_Poisson_1' does it
//         IDBC.clear();
//         gh.clear();
//         tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);

//         // (G) Preserve newly interpolated phi for the next pass
//         //     so we can re-check or re-interpolate from 'phi_new'
//         phi_old = phi_new;
//     }

//     std::cout << "\nAll " << level << " local refinement passes complete.\n";

//     // (H) Optionally output final mesh or other data
//     tt3.OutputGeom_All(path_in + "final_geom");

//     return 0;
// }

int kernel::FindNearestNeighbor(
    const std::vector<BezierElement3D>& bzmesh_old, 
    const BezierElement3D& target_element, 
    double& min_distance, 
    const std::vector<double>& phi_old)
{
    // Initialize the nearest neighbor index and minimum distance
    int nearest_index = -1;  // -1 indicates no valid neighbor found
    min_distance = std::numeric_limits<double>::max();
    
    // Ensure target_element.pts[0] is valid
    if (target_element.pts.empty() || target_element.pts[0].size() < 3) {
        throw std::runtime_error("Invalid target_element: missing point data.");
    }

    // Loop through the old mesh elements to find the nearest neighbor
    for (size_t i = 0; i < bzmesh_old.size(); ++i) {
        // Optionally skip elements with phi_old[i] == 0
        if (phi_old[i] == 0.0) continue; 

        // Ensure the old element has valid point data
        if (bzmesh_old[i].pts.empty() || bzmesh_old[i].pts[0].size() < 3) {
            continue;  // Skip invalid elements
        }

        // Calculate the squared Euclidean distance
        double distance_squared = 0.0;
        for (int j = 0; j < 3; ++j) {
            double diff = target_element.pts[0][j] - bzmesh_old[i].pts[0][j];
            distance_squared += diff * diff;
        }

        // Update the nearest neighbor if a closer one is found
        if (distance_squared < min_distance) {
            min_distance = distance_squared;
            nearest_index = static_cast<int>(i);
        }
    }

    // Convert min_distance to the actual distance
    if (nearest_index != -1) {
        min_distance = std::sqrt(min_distance);
    } else {
        min_distance = 0.0; // No valid neighbors found
    }

    return nearest_index;
}
// // Function to find the index of the nearest neighbor in the old mesh
// int kernel::FindNearestNeighbor(const std::vector<BezierElement3D>& bzmesh_old, const BezierElement3D& target_element, double& min_distance, const std::vector<double>& phi_old)
// {
// 	int nearest_index = 0;
// 	min_distance = std::numeric_limits<double>::max();
// 	// double min_distance = std::numeric_limits<double>::max();

// 	for (int i = 0; i < bzmesh_old.size(); ++i) {
// 		// if (phi_old[i] != 0) {
// 			// Calculate the Euclidean distance between target_element and elements in bzmesh_old
// 			double distance = 0.0;
// 			for (int j = 0; j < 3; ++j)
// 			{
// 				double diff = target_element.pts[0][j] - bzmesh_old[i].pts[0][j];
// 				distance += diff * diff;
// 			}
// 			distance = std::sqrt(distance);

// 			// Update nearest neighbor if a closer one is found
// 			if (distance < min_distance) {
// 				min_distance = distance;
// 				nearest_index = i;
// 			}
// 		// }
// 	}

// 	return nearest_index;
// }

std::vector<double> kernel::InterpolateValues(
    const std::vector<BezierElement3D>& bzmesh_old,
    const std::vector<double>& phi_old,
    const std::vector<BezierElement3D>& bzmesh_new)
{
    std::vector<double> phi_new(bzmesh_new.size(), 0.0);

    for (size_t i = 0; i < bzmesh_new.size(); ++i) {
        double min_distance = std::numeric_limits<double>::max();
        int nearest_index = FindNearestNeighbor(bzmesh_old, bzmesh_new[i], min_distance, phi_old);

        // If an exact match is found, no need for further checks
        if (min_distance <= 0.01) {
            phi_new[i] = phi_old[nearest_index];
            continue;
        }

        // Assign the nearest neighbor value if within the threshold
        const double DIST_THRESHOLD = 4.0; // Threshold for valid interpolation
        if (min_distance <= DIST_THRESHOLD) {
            phi_new[i] = phi_old[nearest_index];
        } else {
            phi_new[i] = 0.0; // Default value for no valid neighbor
        }
    }

    return phi_new;
}
// // Function to perform interpolation from old mesh to new mesh
// std::vector<double> kernel::InterpolateValues(const std::vector<BezierElement3D>& bzmesh_old,
//                                      const std::vector<double>& phi_old,
//                                      const std::vector<BezierElement3D>& bzmesh_new)
// {
// 	std::vector<double> phi_new(bzmesh_new.size(), 0.0);

// 	for (int i = 0; i < bzmesh_new.size(); ++i) {
// 		// Find the nearest neighbor in the old mesh for each element in bzmesh_new
// 		double dist(10);
// 		int nearest_index = FindNearestNeighbor(bzmesh_old, bzmesh_new[i], dist, phi_old);
// 		// std::cout << dist << std::endl;
// 		// Interpolate the value based on the nearest neighbor
// 		// phi_new[i] = phi_old[nearest_index];
// 		if (dist <= 4) {
// 			phi_new[i] = phi_old[nearest_index];
// 		} else {
// 			phi_new[i] = 0;
// 		}
			
// 	}

// 	return phi_new;
// }

void kernel::writeVectorToFile(const std::vector<double>& data, const std::string& filename, bool binary) {
	std::ofstream outfile;

	if (binary) {
		outfile.open(filename, std::ios::out | std::ios::binary);
	} else {
		outfile.open(filename);
	}

	if (!outfile) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}

	if (binary) {
		outfile.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(double));
	} else {
		for (const auto& value : data) {
			outfile << value << " ";
		}
	}

	std::cout << "Vector successfully written to " << filename << std::endl;
	outfile.close();
}

std::vector<double> kernel::readVectorFromFile(const std::string& filename, bool binary) {
	std::ifstream infile;

	if (binary) {
		infile.open(filename, std::ios::in | std::ios::binary);
	} else {
		infile.open(filename);
	}

	if (!infile) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return {};
	}

	std::vector<double> data;

	if (binary) {
		infile.seekg(0, std::ios::end);
		size_t fileSize = infile.tellg();
		infile.seekg(0, std::ios::beg);

		data.resize(fileSize / sizeof(double));
		infile.read(reinterpret_cast<char*>(data.data()), fileSize);
	} else {
		double value;

		while (infile >> value) {
			data.push_back(value);
		}
	}

	std::cout << "Vector successfully read from " << filename << std::endl;
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