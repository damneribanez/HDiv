
#include "path.h"
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZRefPatternTools.h"
#include "CreateDFN.h"

#include "tpzhierarquicalgrid.h"
#include "TPZReadGIDGrid.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "pzgeotetrahedra.h"
#include "pzgeoelrefless.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticline.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"
#include "tpzquadraticquad.h"
#include "tpzarc3d.h"

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "pzfunction.h"
#include "tpzchangeel.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"

#include "TPZPrimalPoisson.h"
#include "pzmatmixedpoisson3d.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZCompMeshTools.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"

#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckmesh.h"
#include "TPZGmshReader.h"
#include "pzintel.h"

#include "TPZHybridizeHDiv.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZLagrangeMultiplier.h"

#include "TFracture.h"
#include "THybridizeDFN.h"
#include "pzl2projection.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "TPZMixedDarcyFlow.h"
#include "TPZTracerFlow.h"
#include "TPZNormalDarcyFlow.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

using namespace std;

struct SimulationCase {
    bool            IsMHMQ;
    bool            IsHybrid;
    bool            UsePardisoQ;
    bool            UseFrontalQ;
    bool            UseGmshMeshQ;
    bool            NonAffineQ;
    int             elemen_type;
    int             n_h_levels;
    int             n_p_levels;
    int             n_acc_terms;
    int             int_order;
    int             n_threads;
    int             perturbation_type;
    std::string     mesh_type;
    std::string     domain_type;
    std::string     conv_summary;
    std::string     dump_folder;
    TPZStack<int>   omega_ids;
    TPZStack<int>   omega_dim;
    TPZStack<int>   gamma_ids;
    TPZStack<int>   gamma_dim;
    TPZStack<REAL>   permeabilities;
    TPZStack<REAL>   porosities;
    TPZStack<REAL>   type;
    TPZStack<REAL>   vals;
    REAL            c_inlet;

    SimulationCase() : IsMHMQ(false), UsePardisoQ(true), IsHybrid(false),UseFrontalQ(false), UseGmshMeshQ(false), NonAffineQ(false), elemen_type(0), n_h_levels(0), n_p_levels(1), n_acc_terms(0), int_order(1), n_threads(0),perturbation_type(0), mesh_type(""), domain_type(""),conv_summary(""),dump_folder(""),omega_ids(),omega_dim(),gamma_ids(), gamma_dim(), permeabilities(),porosities(), type(), vals(), c_inlet(0)
    {

    }

    SimulationCase(const SimulationCase &copy) :IsHybrid(copy.IsHybrid) ,IsMHMQ(copy.IsMHMQ), UsePardisoQ(copy.UsePardisoQ), UseFrontalQ(copy.UseFrontalQ),
    UseGmshMeshQ(copy.UseGmshMeshQ), NonAffineQ(copy.NonAffineQ), elemen_type(copy.elemen_type), n_h_levels(copy.n_h_levels), n_p_levels(copy.n_p_levels), n_acc_terms(copy.n_acc_terms), int_order(copy.int_order),n_threads(copy.n_threads),perturbation_type(copy.perturbation_type), mesh_type(copy.mesh_type), domain_type(copy.domain_type), conv_summary(copy.conv_summary),
    dump_folder(copy.dump_folder), omega_ids(copy.omega_ids), omega_dim(copy.omega_dim), gamma_ids(copy.gamma_ids), gamma_dim(copy.gamma_dim),
    permeabilities(copy.permeabilities),porosities(copy.porosities),
    type(copy.type),
    vals(copy.vals),c_inlet(copy.c_inlet)
    {

    }

    SimulationCase &operator=(const SimulationCase &copy)
    {
        /// check for self-assignment
        if(&copy == this){
            return *this;
        }

        IsMHMQ = copy.IsMHMQ;
        IsHybrid = copy.IsHybrid;
        UsePardisoQ = copy.UsePardisoQ;
        UseFrontalQ = copy.UseFrontalQ;
        UseGmshMeshQ = copy.UseGmshMeshQ;
        NonAffineQ = copy.NonAffineQ;
        elemen_type = copy.elemen_type;
        n_h_levels = copy.n_h_levels;
        n_p_levels = copy.n_p_levels;
        n_acc_terms = copy.n_acc_terms;
        int_order = copy.int_order;
        n_threads = copy.n_threads;
        perturbation_type = copy.perturbation_type;
        mesh_type = copy.mesh_type;
        domain_type = copy.domain_type;
        conv_summary = copy.conv_summary;
        dump_folder = copy.dump_folder;
        omega_ids = copy.omega_ids;
        omega_dim = copy.omega_dim;
        gamma_ids = copy.gamma_ids;
        gamma_dim = copy.gamma_dim;
        permeabilities=copy.permeabilities;
        porosities = copy.porosities;
        type=copy.type;
        vals=copy.vals;
        c_inlet=copy.c_inlet;
        return *this;
    }
};
TPZCompMesh * FluxMesh(TPZGeoMesh * geometry, int order, SimulationCase sim_data);
TPZCompMesh * PressureMesh(TPZGeoMesh * geometry, int order, SimulationCase sim_data);
TPZMultiphysicsCompMesh * MPCMeshMixed(TPZGeoMesh * geometry, int order, SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec);
TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data);
TPZMultiphysicsCompMesh * MPTransportMesh(TPZMultiphysicsCompMesh * mixed, TPZStack<TFracture> & fracture_data ,SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec);
TPZCompMesh * SaturationMesh(TPZMultiphysicsCompMesh *mpcompmesh, int order, SimulationCase sim);
void InsertTransportInterfaceElements(TPZMultiphysicsCompMesh *cmesh);
void InsertInterfacesBetweenElements(int transport_matid, TPZCompMesh * cmesh, std::vector<int> & cel_indexes);
TPZCompMesh *CreateTransportMesh(TPZMultiphysicsCompMesh *cmesh, TPZStack<TFracture> &fractures);
void CreateTransportElement(int p_order, TPZCompMesh *cmesh, TPZGeoEl *gel, bool is_BC);
TPZAnalysis * CreateTransportAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data);
TPZFMatrix<STATE> TimeForward(TPZAnalysis * tracer_analysis, int & n_steps, REAL & dt, TPZFMatrix<STATE> & M_diag);
void UniformRefinement(TPZGeoMesh * geometry, int h_level) ;
/// Executes cube
void Pretty_cube();
void Pretty_cubeII();
void Case_1();
void TesisDamner();

int main(){
    
    std::string source_dir = SOURCE_DIR;
#ifdef LOG4CXX
    std::string log_file = source_dir;
    log_file += "/dfn.cfg";
    InitializePZLOG(log_file);
#endif
    
     TesisDamner();
 //   Pretty_cubeII();
}
void Pretty_cube(){
    
    
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.IsHybrid=true;
    sim.n_threads = 8;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0);
    sim.porosities.push_back(1.0);
    
    /// not used but inserted
    sim.omega_ids.push_back(2);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0);
    sim.porosities.push_back(1.0);
    
    /// C inlet value
    sim.c_inlet = 1.0;
    
    int bc_inlet  = 3;
    int bc_outlet = 4;
    int bc_non_flux = 5;
    
    sim.gamma_ids.push_back(bc_inlet);
    sim.gamma_dim.push_back(3);
    sim.gamma_ids.push_back(bc_outlet);
    sim.gamma_dim.push_back(3);
    sim.gamma_ids.push_back(bc_non_flux);
    sim.gamma_dim.push_back(3);
    
    int bc_type_D = 0;    //    D = 0;
    int bc_type_N = 1;    //    N = 1;
    REAL p_inlet  = 2.0;
    REAL p_outlet = 1.0;
    REAL qn       = 0.0;
    
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_N);
    
    sim.vals.push_back(p_inlet);
    sim.vals.push_back(p_outlet);
    sim.vals.push_back(qn);
    
    /// Defining DFN boundary data (id,bc_type,data)
    std::vector<std::tuple<int,int,REAL>> bc_ids_2d;
    bc_ids_2d.push_back(std::make_tuple(bc_inlet,bc_type_D,p_inlet));
    bc_ids_2d.push_back(std::make_tuple(bc_outlet,bc_type_D,p_outlet));
    bc_ids_2d.push_back(std::make_tuple(bc_non_flux,bc_type_N,qn));
    
    int bc_1d_inlet  = 130;
    int bc_1d_outlet = 140;
    int bc_1d_non_flux = 150;
    
    int bc_0d_inlet  = 230;
    int bc_0d_outlet = 240;
    int bc_0d_non_flux = 250;
    
    std::map<int,int> bc_ids_1d_map;
    bc_ids_1d_map.insert(std::make_pair(bc_inlet,bc_1d_inlet));
    bc_ids_1d_map.insert(std::make_pair(bc_outlet,bc_1d_outlet));
    bc_ids_1d_map.insert(std::make_pair(bc_non_flux,bc_1d_non_flux));
    
    std::map<int,int> bc_ids_0d_map;
    bc_ids_0d_map.insert(std::make_pair(bc_inlet,bc_0d_inlet));
    bc_ids_0d_map.insert(std::make_pair(bc_outlet,bc_0d_outlet));
    bc_ids_0d_map.insert(std::make_pair(bc_non_flux,bc_0d_non_flux));
    
    /// Defining DFN data
    TPZStack<TFracture> fracture_data;
    TFracture fracture;
    fracture.m_id.insert(6);
    fracture.m_id.insert(7);
    fracture.m_id.insert(8);
    fracture.m_dim              = 2;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 0.5;
    fracture_data.push_back(fracture);
    fracture.m_id.clear();
    fracture.m_id.insert(9);
    fracture.m_dim              = 1;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 1.0;
    fracture_data.push_back(fracture);
    fracture.m_id.clear();
    fracture.m_id.insert(10);
    fracture.m_dim              = 0;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 1.0;
    fracture_data.push_back(fracture);
    
    
    /// Benchmarks Material ID convention
    /// 1 and 2 for 3D matrix
    /// 3,4, and 5 for 2D matrix boundaries (3 -> inlet, 4 -> outlet, 5 -> impervious)
    /// 6 fractures
    /// 7 fractures intersections
    /// 8 crossing intersections
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag(4); // From 0D to 3D
    dim_name_and_physical_tag[3]["RockMatrix_1"] = 1;
    dim_name_and_physical_tag[3]["RockMatrix_2"] = 2;
    dim_name_and_physical_tag[2]["BCInlet"] = 3;
    dim_name_and_physical_tag[2]["BCOutlet"] = 4;
    dim_name_and_physical_tag[2]["BCImpervious"] = 5;
    dim_name_and_physical_tag[2]["Fractures_1"] = 6;
    dim_name_and_physical_tag[2]["Fractures_2"] = 7;
    dim_name_and_physical_tag[2]["Fractures_3"] = 8;
    dim_name_and_physical_tag[1]["FracturesIntersections"] = 9;
    dim_name_and_physical_tag[0]["CrossingIntresections"] = 10;
    
    TPZGmshReader Geometry;
    std::string source_dir = SOURCE_DIR;
    std::string file_gmsh = source_dir + "/meshes/the_cuttest_cube/cube.msh";
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    std::string version("4.1");
    Geometry.SetFormatVersion(version);
    
    Geometry.SetDimNamePhysical(dim_name_and_physical_tag);
    gmesh = Geometry.GeometricGmshMesh(file_gmsh.c_str());
    Geometry.PrintPartitionSummary(std::cout);
    
    
#ifdef PZDEBUG
    std::ofstream file("geometry_cube_base.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    std::ofstream file_txt("geometry_cube_base.txt");
    gmesh->Print(file_txt);
#endif
    
    int p_order = 1;
    TPZVec<TPZCompMesh *> meshvec;
    TPZCompMesh *cmixedmesh = NULL;
    cmixedmesh = MPCMeshMixed(gmesh, p_order, sim, meshvec);
#ifdef PZDEBUG
    std::ofstream filemixed("mixed_cmesh.txt");
  //  cmixedmesh->Print(filemixed);
#endif
    
    TPZCompMesh *cmeshm =NULL;
    THybridizeDFN dfn_hybridzer;
    dfn_hybridzer.SetFractureData(fracture_data);
    dfn_hybridzer.SetReservoirBoundaryData(bc_ids_2d);
    dfn_hybridzer.SetMapReservoirBCToDFNBC1DIds(bc_ids_1d_map);
    dfn_hybridzer.SetMapReservoirBCToDFNBC0DIds(bc_ids_0d_map);
    cmeshm = dfn_hybridzer.Hybridize(cmixedmesh);
    
    
    TPZMultiphysicsCompMesh * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmeshm);
    TPZManVector<TPZCompMesh *> vectormesh = mp_cmesh->MeshVector();
    TPZCompMesh *Flux_Mesh = vectormesh[0] ;
    TPZCompMesh *Pressure_Mesh = vectormesh[1] ;
    std::ofstream file3("PRESSURE.vtk");
     std::ofstream file2("FLUX.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(Flux_Mesh, file2);
    TPZVTKGeoMesh::PrintCMeshVTK(Pressure_Mesh, file3);
    
#ifdef PZDEBUG
    {
        std::ofstream file("geometry_cube_after.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
        std::ofstream file_txt("geometry_cube_after.txt");
        gmesh->Print(file_txt);
    }
#endif
    
    
   // TPZCompMesh *s_cmesh = CreateTransportMesh(mp_cmesh, fracture_data);
    
#ifdef PZDEBUG
    {
//        std::ofstream file("transport_mesh.vtk");
//        TPZVTKGeoMesh::PrintCMeshVTK(s_cmesh, file);
//        std::ofstream out("transport_mesh.txt");
//        s_cmesh->Print(out);
    }
#endif
    
//    TPZManVector<TPZCompMesh *,3> meshtrvec(3);
//    meshtrvec[0] = meshvec[0];
//    meshtrvec[1] = meshvec[1];
//    meshtrvec[2] = s_cmesh;
    
   // TPZMultiphysicsCompMesh *cmesh_transport = MPTransportMesh(mp_cmesh, fracture_data, sim, meshtrvec);
  //  TPZAnalysis * tracer_analysis = CreateTransportAnalysis(cmesh_transport, sim);
    
    bool solve_dfn_problem_Q = true;
    if (solve_dfn_problem_Q) {
        TPZManVector<TPZCompMesh * > mesh_vec = mp_cmesh->MeshVector();
        
#ifdef PZDEBUG
        {
            std::ofstream file_hybrid_mixed_q("Hybrid_mixed_cmesh_q.txt");
            mesh_vec[0]->ComputeNodElCon();
            mesh_vec[0]->Print(file_hybrid_mixed_q);
            
            std::ofstream file_hybrid_mixed_p("Hybrid_mixed_cmesh_p.txt");
            mesh_vec[1]->ComputeNodElCon();
            mesh_vec[1]->Print(file_hybrid_mixed_p);
            
            std::ofstream file_hybrid_mixed("Hybrid_mixed_cmesh.txt");
            cmeshm->ComputeNodElCon();
            cmeshm->Print(file_hybrid_mixed);
        }
#endif
        
        std::cout << "Condensing DFN equations." << std::endl;
        std::cout << "DFN neq before condensation = " << mp_cmesh->NEquations() << std::endl;
        dfn_hybridzer.GroupElements(mp_cmesh);
        std::cout << "DFN neq = " << mp_cmesh->NEquations() << std::endl;
        
#ifdef PZDEBUG
        {
            
            std::ofstream file_hybrid_condensed("Hybrid_mixed_condensed.txt");
            mp_cmesh->ComputeNodElCon();
            mp_cmesh->Print(file_hybrid_condensed);
        }
#endif
        
        std::cout << "DFN equations are condensed." << std::endl;
        
        TPZAnalysis *an = CreateAnalysis(mp_cmesh, sim);
        std::cout << "Assembly DFN problem neq = " << mp_cmesh->NEquations() << std::endl;
        an->Assemble();
        std::cout << "Assembly for DFN complete." << std::endl;
        
        std::cout << "Solving DFN problem." << std::endl;
        an->Solve();
        std::cout << "DFN problem solved." << std::endl;
        mp_cmesh->LoadSolutionFromMultiPhysics();
        
#ifdef PZDEBUG
        std::ofstream file_geo_hybrid("geometry_cube_hybrid.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(cmeshm->Reference(), file_geo_hybrid);
#endif
        
        TPZStack<std::string,10> scalnames, vecnames;
        vecnames.Push("q");
        scalnames.Push("p");
        
        int div = 0;
        std::set<int> mat_id_3D;
        mat_id_3D.insert(1);
        mat_id_3D.insert(2);
        std::string file_reservoir("cube.vtk");
        an->DefineGraphMesh(3,mat_id_3D,scalnames,vecnames,file_reservoir);
        an->PostProcess(div,3);
        
        std::set<int> mat_id_2D;
        mat_id_2D.insert(6);
        mat_id_2D.insert(7);
        mat_id_2D.insert(8);
        std::string file_frac("fracture.vtk");
        an->DefineGraphMesh(2,mat_id_2D,scalnames,vecnames,file_frac);
        an->PostProcess(div,2);
        
        std::set<int> mat_id_1D;
        mat_id_1D.insert(7);
         mat_id_1D.insert(8);
         mat_id_1D.insert(9);
        std::string file_frac_intersections("fracture_intersections.vtk");
        an->DefineGraphMesh(1,mat_id_1D,scalnames,vecnames,file_frac_intersections);
        an->PostProcess(div,1);
    }
    
    return;
}
TPZMultiphysicsCompMesh * MPCMeshMixed(TPZGeoMesh * geometry, int order, SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec){
    
    int dimension = geometry->Dimension();
    int nvols = sim_data.omega_ids.size();
    int nbound= sim_data.gamma_ids.size();
    if (nvols<1) {
        std::cout<<"Error: Omega is not defined."<<std::endl;
        DebugStop();
    }
    if (nbound<1) {
        std::cout<<"Error: Gamma is not defined."<<std::endl;
        DebugStop();
    }
    
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(geometry);
    TPZFNMatrix<9,STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        
        TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(sim_data.omega_ids[ivol],dimension);
        volume->SetPermeability(sim_data.permeabilities[ivol]);
        cmesh->InsertMaterialObject(volume);
        
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                val2(0,0)=sim_data.vals[ibound];
                int condType=sim_data.type[ibound];
                TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
            }
        }
    }
    
    cmesh->SetDimModel(dimension);
    TPZManVector<TPZCompMesh * ,2> mesh_vec(2);
    /// generate the 3D flux mesh (only)
    mesh_vec[0] = FluxMesh(geometry, order, sim_data);
    mesh_vec[1] = PressureMesh(geometry, order, sim_data);
    TPZManVector<int,5> active_approx_spaces(2); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    cmesh->BuildMultiphysicsSpace(active_approx_spaces,mesh_vec);
    
    
    
    if (sim_data.IsMHMQ) {
        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
    }
    //    else{
    //        TPZCompMeshTools::GroupElements(cmesh);
    //        std::cout << "Created grouped elements\n";
    //        bool keepmatrix = false;
    //        bool keeponelagrangian = true;
    //        TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);
    //        std::cout << "Created condensed elements\n";
    //        cmesh->CleanUpUnconnectedNodes();
    //        cmesh->ExpandSolution();
    //    }
    
    std::cout << "Created multi-physics DFN mesh\n";
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name  << "Dual_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    meshvec = mesh_vec;
    return cmesh;
}
TPZCompMesh * FluxMesh(TPZGeoMesh * geometry, int order, SimulationCase sim_data){
    
    int dimension = geometry->Dimension();
    int nvols = sim_data.omega_ids.size();
    int nbound = sim_data.gamma_ids.size();
    
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(sim_data.omega_ids[ivol],sim_data.omega_dim[ivol]);
        volume->SetPermeability(sim_data.permeabilities[ivol]);
        
        cmesh->InsertMaterialObject(volume);
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                if(sim_data.gamma_dim[ibound] != sim_data.omega_dim[ivol]) continue;
                val2(0,0)=sim_data.vals[ibound];
                int condType=sim_data.type[ibound];
                TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
            }
        }
    }
    
//    for (int dim=3; dim>0; dim--) {
//        cmesh->Reference()->ResetReference();
//        cmesh->SetDimModel(dim);
//        cmesh->SetDefaultOrder(order);
//        std::set<int> matids;
//        for (int ivol=0; ivol<nvols; ivol++) {
//            if(sim_data.omega_dim[ivol]==dim) matids.insert(sim_data.omega_ids[ivol]);
//        }
//        for (int ibound=0; ibound<nbound; ibound++) {
//            if(sim_data.gamma_dim[ibound]==dim) matids.insert(sim_data.gamma_ids[ibound]);
//        }
//        cmesh->SetAllCreateFunctionsHDiv();
//        cmesh->AutoBuild(matids);
//    }
    
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    
    cmesh->SetDimModel(geometry->Dimension());
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name << "q_cmesh_raw" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
    
}
TPZCompMesh * PressureMesh(TPZGeoMesh * geometry, int order, SimulationCase sim_data){
    
    int dimension = geometry->Dimension();
    int nvols = sim_data.omega_ids.size();
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    std::set<int> matids;
    for (int ivol=0; ivol < nvols; ivol++) {
        TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(sim_data.omega_ids[ivol],dimension);
        volume->SetPermeability(sim_data.permeabilities[ivol]);
        cmesh->InsertMaterialObject(volume);
        if (sim_data.omega_dim[ivol] == dimension) {
            matids.insert(sim_data.omega_ids[ivol]);
        }
    }
    
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    cmesh->AutoBuild(matids);
    cmesh->InitializeBlock();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name << "p_cmesh_raw" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
    
}
TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh, true);
    
    if (sim_data.UsePardisoQ) {
        
        TPZSymetricSpStructMatrix matrix(cmesh);
        //        TPZSkylineStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);
        
        return analysis;
    }
    
    if (sim_data.UseFrontalQ) {
        
        TPZParFrontStructMatrix<TPZFrontSym<STATE> > matrix(cmesh);
        matrix.SetDecomposeType(ELDLt);
        matrix.SetNumThreads(sim_data.n_threads);
        
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);
        
        return analysis;
    }
    else{
        
        TPZSkylineStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ECholesky);
        analysis->SetSolver(step);
        analysis->SetStructuralMatrix(matrix);
        return analysis;
    }
    
    return analysis;
    
}
void Pretty_cubeII(){
    
    //Lectura de la malha
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag(4); // From 0D to 3D
    dim_name_and_physical_tag[2]["RockMatrix_1"] = 1;
    dim_name_and_physical_tag[2]["RockMatrix_2"] = 2;
    dim_name_and_physical_tag[1]["BCInlet"] = 3;
    dim_name_and_physical_tag[1]["BCOutlet"] = 4;
    dim_name_and_physical_tag[1]["BCImpervious"] = 5;
    dim_name_and_physical_tag[1]["Fractures_1"] = 6;
    dim_name_and_physical_tag[1]["Fractures_2"] = 6;
    dim_name_and_physical_tag[1]["Fractures_3"] = 6;
    dim_name_and_physical_tag[0]["FracturesIntersections"] = 7;
    dim_name_and_physical_tag[0]["CrossingIntresections"] = 8;
    
    TPZGmshReader Geometry;
    std::string source_dir = SOURCE_DIR;
    std::string file_gmsh = source_dir + "/meshes/the_cuttest_cube/cube.msh";
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    std::string version("4.1");
    Geometry.SetFormatVersion(version);
    
    Geometry.SetDimNamePhysical(dim_name_and_physical_tag);
    gmesh = Geometry.GeometricGmshMesh(file_gmsh.c_str());
    Geometry.PrintPartitionSummary(std::cout);
    
    UniformRefinement(gmesh, 2);
    
#ifdef PZDEBUG
    std::ofstream file("geometry_cube_base.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    std::ofstream file_txt("geometry_cube_base.txt");
    gmesh->Print(file_txt);
    gmesh->SetDimension(3);
#endif
    
    //Configuracion de Materiales
    
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.IsHybrid=true;
    sim.n_threads = 1;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0);
    sim.porosities.push_back(1.0);
    
    /// not used but inserted
    sim.omega_ids.push_back(2);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0);
    sim.porosities.push_back(1.0);
    
    /// C inlet value
    sim.c_inlet = 1.0;
    
    int bc_inlet  = 3;
    int bc_outlet = 4;
    int bc_non_flux = 5;
    
    sim.gamma_ids.push_back(bc_inlet);
    sim.gamma_dim.push_back(3);
    sim.gamma_ids.push_back(bc_outlet);
    sim.gamma_dim.push_back(3);
    sim.gamma_ids.push_back(bc_non_flux);
    sim.gamma_dim.push_back(3);
    
    int bc_type_D = 0;    //    D = 0;
    int bc_type_N = 1;    //    N = 1;
    REAL p_inlet  = 2.0;
    REAL p_outlet = 1.5;
    REAL qn       = 0.0;
    
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_N);
    
    sim.vals.push_back(p_inlet);
    sim.vals.push_back(p_outlet);
    sim.vals.push_back(qn);
    TPZVec<TPZCompMesh *> meshvec;
    
    /// Defining DFN data
    TPZStack<TFracture> fracture_data;
    TFracture fracture;
    fracture.m_id.insert(6);
    fracture.m_dim              = 2;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 0.5;
    fracture_data.push_back(fracture);
    fracture.m_id.clear();
    fracture.m_id.insert(7);
    fracture.m_dim              = 1;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 1.0;
    fracture_data.push_back(fracture);
    fracture.m_id.clear();
    fracture.m_id.insert(8);
    fracture.m_dim              = 0;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 1.0;
    fracture_data.push_back(fracture);
 
    /// Defining DFN boundary data (id,bc_type,data)
    std::vector<std::tuple<int,int,REAL>> bc_ids_2d;
    bc_ids_2d.push_back(std::make_tuple(bc_inlet,bc_type_D,p_inlet));
    bc_ids_2d.push_back(std::make_tuple(bc_outlet,bc_type_D,p_outlet));
    bc_ids_2d.push_back(std::make_tuple(bc_non_flux,bc_type_N,qn));
    
    int bc_1d_inlet  = 130;
    int bc_1d_outlet = 140;
    int bc_1d_non_flux = 150;
    
    int bc_0d_inlet  = 230;
    int bc_0d_outlet = 240;
    int bc_0d_non_flux = 250;
    
    std::map<int,int> bc_ids_1d_map;
    bc_ids_1d_map.insert(std::make_pair(bc_inlet,bc_1d_inlet));
    bc_ids_1d_map.insert(std::make_pair(bc_outlet,bc_1d_outlet));
    bc_ids_1d_map.insert(std::make_pair(bc_non_flux,bc_1d_non_flux));
    
    std::map<int,int> bc_ids_0d_map;
    bc_ids_0d_map.insert(std::make_pair(bc_inlet,bc_0d_inlet));
    bc_ids_0d_map.insert(std::make_pair(bc_outlet,bc_0d_outlet));
    bc_ids_0d_map.insert(std::make_pair(bc_non_flux,bc_0d_non_flux));
    //
    CreateDFN dfn;
    int flux_trace_id=10;
    int lagrange_id=11;
    int flux_resistivity_id=13;
    int mp_nterface_id =12;
    
    dfn.SetFractureData(fracture_data);
    dfn.SetReservoirBoundaryData(bc_ids_2d);
    dfn.SetMapReservoirBCToDFNBC1DIds(bc_ids_1d_map);
    dfn.SetMapReservoirBCToDFNBC0DIds(bc_ids_0d_map);
    
  //  UniformRefinement(gmesh, 1);
    TPZMultiphysicsCompMesh *mp_initial_cmesh = MPCMeshMixed(gmesh, 1, sim, meshvec);
    dfn.SetPeriferalMaterialIds(flux_trace_id, lagrange_id, mp_nterface_id,flux_resistivity_id);
    TPZCompMesh * dfn_hybrid_cmesh = dfn.CreateDFNCmesh(mp_initial_cmesh);
    
    TPZMultiphysicsCompMesh * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(dfn_hybrid_cmesh);
    
    
    TPZCompMesh *s_cmesh = CreateTransportMesh(mp_cmesh, fracture_data);
    std::ofstream printsat("s_cmesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(s_cmesh, printsat);
    
    std::ofstream printsat2("s_cmesh.txt");
    s_cmesh->Print(printsat2);
    
    TPZManVector<TPZCompMesh *,3> meshtrvec(3);
    meshtrvec[0] = meshvec[0];
    meshtrvec[1] = meshvec[1];
    meshtrvec[2] = s_cmesh;
    
    TPZMultiphysicsCompMesh *cmesh_transport = MPTransportMesh(mp_cmesh, fracture_data, sim, meshtrvec);
    TPZAnalysis * tracer_analysis = CreateTransportAnalysis(cmesh_transport, sim);

    
    bool solve_dfn_problem_Q = true;
    if (solve_dfn_problem_Q) {
        TPZManVector<TPZCompMesh * > mesh_vec = mp_cmesh->MeshVector();
        
        std::ofstream printpress("pres.vtk");
  //      mesh_vec[1]->Print(printpress);
       TPZVTKGeoMesh::PrintCMeshVTK(mesh_vec[1], printpress);

        std::ofstream printflux("flux.vtk");
  //      mesh_vec[0]->Print(printflux);
       TPZVTKGeoMesh::PrintCMeshVTK(mesh_vec[0], printflux);
        
        std::cout << "Condensing DFN equations." << std::endl;
        std::cout << "DFN neq before condensation = " << mp_cmesh->NEquations() << std::endl;
      //  dfn.GroupElements(mp_cmesh);
        std::cout << "DFN neq = " << mp_cmesh->NEquations() << std::endl;
        
#ifdef PZDEBUG
        {
            
//            std::ofstream file_hybrid_condensed("Hybrid_mixed_condensed.txt");
//            mp_cmesh->ComputeNodElCon();
//            mp_cmesh->Print(file_hybrid_condensed);
        }
#endif
        
        std::cout << "DFN equations are condensed." << std::endl;
        
        TPZAnalysis *an = CreateAnalysis(mp_cmesh, sim);
        std::cout << "Assembly DFN problem neq = " << mp_cmesh->NEquations() << std::endl;
        an->Assemble();
        std::cout << "Assembly for DFN complete." << std::endl;
        
        std::cout << "Solving DFN problem." << std::endl;
        an->Solve();
        std::cout << "DFN problem solved." << std::endl;
        mp_cmesh->LoadSolutionFromMultiPhysics();
        
        
        TPZStack<std::string,10> scalnames, vecnames;
        vecnames.Push("q");
        scalnames.Push("p");
        
        int div = 0;
        std::set<int> mat_id_3D;
        mat_id_3D.insert(1);
//        mat_id_3D.insert(2);
        std::string file_reservoir("cube.vtk");
        an->DefineGraphMesh(3,mat_id_3D,scalnames,vecnames,file_reservoir);
        an->PostProcess(div,3);
        
        std::set<int> mat_id_2D;
        mat_id_2D.insert(6);
        std::string file_frac("fracture.vtk");
        an->DefineGraphMesh(2,mat_id_2D,scalnames,vecnames,file_frac);
        an->PostProcess(div,2);
        
        std::set<int> mat_id_1D;
        mat_id_1D.insert(7);
        std::string file_frac_intersections("fracture_intersections.vtk");
        an->DefineGraphMesh(1,mat_id_1D,scalnames,vecnames,file_frac_intersections);
        an->PostProcess(div,1);
    }
    
    int n_steps = 10;
    REAL dt     = 10.0;
    TPZFMatrix<STATE> M_diag;
    TPZFMatrix<STATE> saturations = TimeForward(tracer_analysis, n_steps, dt, M_diag);
    
    
    
    
    return;
}
TPZMultiphysicsCompMesh * MPTransportMesh(TPZMultiphysicsCompMesh * mixed, TPZStack<TFracture> & fracture_data ,SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec){
    
    
    TPZGeoMesh *geometry = mixed->Reference();
    int dimension = geometry->Dimension();
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(geometry);
    
    /// Inserting matrix materials
    int n_vols = sim_data.omega_ids.size();
    for (int i = 0; i < n_vols; i++) {
        int mat_id = sim_data.omega_ids[i];
        REAL phi = sim_data.porosities[i];
        TPZTracerFlow * volume = new TPZTracerFlow(mat_id,0);
        volume->SetPorosity(phi);
        cmesh->InsertMaterialObject(volume);
    }
    
    /// Inserting fracture materials
    int n_fracs = fracture_data.size();
    for (int i = 0; i < n_fracs; i++) {
        for(auto mat_id :fracture_data[i].m_id)
        {
            REAL phi = fracture_data[i].m_porosity;
            REAL d_opening = fracture_data[i].m_d_opening;
            TPZTracerFlow * volume = new TPZTracerFlow(mat_id,0);
            volume->SetPorosity(phi);
            volume->SetFractureCrossLength(d_opening);
            cmesh->InsertMaterialObject(volume);
        }
    }
    
    std::set<int> bc_inlet_mat_ids = {3,5,130,150,230,250,-1942};
    std::set<int> bc_outlet_mat_ids = {4,140,240};
    
    TPZMaterial * material = cmesh->FindMaterial(1);
    if (!material) {
        DebugStop();
    }
    TPZFMatrix<STATE> val1(1,1,0);
    TPZFMatrix<STATE> val2(1,1,0);
    
    int typ_inlet = 0; // inlet
    /// Inserting the materials
    val2(0,0) = sim_data.c_inlet;
    for (auto mat_id: bc_inlet_mat_ids) {
        TPZMaterial * bc = material->CreateBC(material, mat_id, typ_inlet, val1, val2);
        cmesh->InsertMaterialObject(bc);
    }
    
    int typ_outlet = 1; // outlet
    /// Inserting the materials
    val2(0,0) = 0.0;
    for (auto mat_id: bc_outlet_mat_ids) {
        TPZMaterial * bc = material->CreateBC(material, mat_id, typ_outlet, val1, val2);
        cmesh->InsertMaterialObject(bc);
    }
    
    cmesh->SetDimModel(dimension);
    
    TPZManVector<int,5> active_approx_spaces(3); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 0;
    active_approx_spaces[1] = 0;
    active_approx_spaces[2] = 1;
    cmesh->BuildMultiphysicsSpace(active_approx_spaces,meshvec);
    
    
    if (sim_data.IsMHMQ) {
        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
    }
    else{
        //        TPZCompMeshTools::GroupElements(cmesh);
        //        std::cout << "Created grouped elements\n";
        //        bool keepmatrix = false;
        //        bool keeponelagrangian = true;
        //        TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);
        //        std::cout << "Created condensed elements\n";
        //        cmesh->CleanUpUnconnectedNodes();
        //        cmesh->ExpandSolution();
    }
    
    std::ofstream filesat("saturation.vtk");
   TPZVTKGeoMesh::PrintGMeshVTK(geometry, filesat);
   InsertTransportInterfaceElements(cmesh);
    
    std::cout << "Created multi-physics transport mesh\n";
    
#ifdef PZDEBUG
    std::ofstream transport("transport_cmesh.txt");
    cmesh->Print(transport);
#endif
    
    return cmesh;
}
TPZCompMesh * SaturationMesh(TPZMultiphysicsCompMesh *mpcompmesh, int order, SimulationCase sim){
    
    TPZManVector<TPZCompMesh * > mesh_vec = mpcompmesh->MeshVector();
    TPZCompMesh *qMesh = mesh_vec[0];
    TPZCompMesh *pMesh = mesh_vec[1];
    mpcompmesh->MeshVector().Resize(3);
    
    int dim =3;
    pMesh->LoadReferences();
    
    TPZGeoMesh *gmesh = pMesh->Reference();
    TPZCompMesh *smesh = new TPZCompMesh(gmesh);
    int mat_id_s = 1;
    
    
    int nel = pMesh->NElements();
    TPZStack<TPZCompElSide> el_oned_connec;
    for(int iel =0; iel<nel; iel++){
        TPZGeoEl *gel = pMesh->Element(iel)->Reference();
        if (gel->Dimension()==0) {
            TPZGeoElSide gelside(gel,0);
            TPZStack<TPZCompElSide> elsidevec;
            int onlyinterpolated =0;
            int removeduplicates =0;
            gelside.EqualLevelCompElementList(elsidevec, onlyinterpolated,  removeduplicates);
            int nelconnect = elsidevec.size();
            for (int ioned=0; ioned<nelconnect; ioned++) {
                if (nelconnect>2) {
                    int64_t index;
                    gel->SetMaterialId(mat_id_s);
                    smesh->CreateCompEl(gel, index);
                }
            }
        }
    }
    
    qMesh->LoadReferences();
    qMesh->SetDimModel(dim);
    
    int nelq = qMesh->NElements();
    for (int Dims=3; Dims<1; Dims++) {
        for (int iel =0; iel<nelq; iel++) {
            TPZGeoEl *gel = qMesh->Element(iel)->Reference();
            if (!gel->Reference()) {
                continue;
            }
            int dimel = gel->Dimension();
            int nconnecs = gel->Reference()->NConnects();
            if (dimel == Dims && nconnecs ==1) {
                int64_t index;
                gel->SetMaterialId(mat_id_s);
                smesh->CreateCompEl(gel, index);
                
            }
        }
    }
    
    
    smesh->SetDimModel(dim);
    smesh->SetDefaultOrder(order);
    
    smesh->SetAllCreateFunctionsDiscontinuous();
    //smesh->AutoBuild();
    
    smesh->AdjustBoundaryElements();
    smesh->CleanUpUnconnectedNodes();
    
    mpcompmesh->MeshVector()[2]=smesh;
    
    return smesh;
    
}
void InsertTransportInterfaceElements(TPZMultiphysicsCompMesh *cmesh)
{
    
    int transport_matid = 1;
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    int64_t nel = cmesh->NElements();
    TPZManVector<std::vector<int>,4> gel_vol_index_to_cel(4);
    
    for (int64_t el = 0; el< nel; el++) {
        
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) DebugStop();
        TPZMultiphysicsElement *celmp = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if(!celmp) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        
        int gel_dim = gel->Dimension();
        gel_vol_index_to_cel[gel_dim].push_back(el);
        
    }
 
    InsertInterfacesBetweenElements(transport_matid, cmesh, gel_vol_index_to_cel[3]);
    InsertInterfacesBetweenElements(transport_matid, cmesh, gel_vol_index_to_cel[2]);
    //InsertInterfacesBetweenElements(transport_matid, cmesh, gel_vol_index_to_cel[1]);
    
    cmesh->ComputeNodElCon();
}

void InsertInterfacesBetweenElements(int transport_matid, TPZCompMesh * cmesh, std::vector<int> & cel_indexes){
    
    TPZGeoMesh * geometry = cmesh->Reference();
    if (!geometry) {
        DebugStop();
    }
    
    int mesh_dim = geometry->Dimension();
    
    std::set<int> bcmat_ids = {3,4,5,130,140,150,230,240,250,-1942};
    bool needs_all_boundaries_Q = false;
    TPZManVector<int64_t,3> left_mesh_indexes(2,0);
    left_mesh_indexes[0] = 0;
    left_mesh_indexes[1] = 2;
    TPZManVector<int64_t,3> right_mesh_indexes(1,0);
    right_mesh_indexes[0] = 2;
    for (auto cel_index: cel_indexes) {
        TPZCompEl *cel = cmesh->Element(cel_index);
        TPZGeoEl *gel = cel->Reference();
        
        int nsides = gel->NSides();
        int geldim = gel->Dimension();
        for (int side = 0; side<nsides; side++) {
            int sidedim = gel->SideDimension(side);
            if(sidedim < geldim-1) continue;
            int mat_id = gel->MaterialId();
            if (bcmat_ids.find(mat_id) != bcmat_ids.end()) {
                continue;
            }
            TPZStack<TPZCompElSide> celstack;
            TPZCompElSide celside(cel,side);
            TPZGeoElSide gelside(gel,side);
            gelside.EqualLevelCompElementList(celstack, 0, 0);
            if(celstack.size() == 0 && sidedim < mesh_dim && needs_all_boundaries_Q)
            {
                DebugStop();
            }
            if(celstack.size() == 0){
                continue;
            }
            
            if(sidedim == geldim-1)
            {
                
                int count_dim_m_1 = 0;
                int stack_index_dim_m_1 = 0;
                int count = 0;
                int stack_index = 0;
                bool same_dimension_Q = false;
                for (int icel = 0; icel < celstack.size(); icel++) {
                    TPZCompEl * cel_neigh = celstack[icel].Element();
                    TPZMultiphysicsInterfaceElement * mp_interface_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel_neigh);
                    int gel_neigh_dim = cel_neigh->Dimension();
                    if (!mp_interface_cel) {
                        if (gel_neigh_dim == geldim-1) {
                            count_dim_m_1++;
                            stack_index_dim_m_1 = icel;
                        }else if (gel_neigh_dim == geldim){
                            count++;
                            stack_index = icel;
                            same_dimension_Q = true;
                        }
                    }
                }
                
                if (same_dimension_Q && count_dim_m_1 == 0) {
                    
                    // There must be a lower dimensional transport element to equate the fluxes
                    if(count != 1){
                        DebugStop();
                    }
                    
                    TPZGeoEl *neighgel = celstack[stack_index].Element()->Reference();
                    if(geldim < neighgel->Dimension())
                    {
                        // we only create interfaces from lower to higher dimensional elements
                        // neighgel must be a boundary element
                        continue;
                    }
                    if(geldim == neighgel->Dimension() && gel->Id() > neighgel->Id())
                    {
                        // the interface is created from the lower id to the higher id
                        continue;
                    }
                    
                    // we need to create the interface
                    TPZGeoElBC gbc(gelside,transport_matid);
                    int64_t index;
                    TPZMultiphysicsInterfaceElement * mp_interface_el = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside, celstack[stack_index]);
                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                    continue;
                }
                
                /// There must be a lower dimensional transport element to equate the fluxes
                if(count_dim_m_1 != 1){
                    DebugStop();
                }
//
                TPZGeoEl *neighgel = celstack[stack_index_dim_m_1].Element()->Reference();
                if(geldim < neighgel->Dimension())
                {
                    // we only create interfaces from lower to higher dimensional elements
                    // neighgel must be a boundary element
                    continue;
                }
                if(geldim == neighgel->Dimension() && gel->Id() > neighgel->Id())
                {
                    // the interface is created from the lower id to the higher id
                    continue;
                }
                {
                    int mat_interface_id;
                    int neigh_mat_id = neighgel->MaterialId();
                    if (bcmat_ids.find(neigh_mat_id) != bcmat_ids.end()) {
                        mat_interface_id = neigh_mat_id;
                    }else{
                        mat_interface_id = transport_matid;
                    }
                    
                    // we need to create the interface
                    TPZGeoElBC gbc(gelside,mat_interface_id);
                    int64_t index;
                    TPZMultiphysicsInterfaceElement * mp_interface_el = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside, celstack[stack_index_dim_m_1]);
                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                }
            }
        }
    }
    
}


TPZCompMesh *CreateTransportMesh(TPZMultiphysicsCompMesh *cmesh, TPZStack<TFracture> &fractures)
{
    TPZCompMesh *q_cmesh = cmesh->MeshVector()[0];
    TPZGeoMesh * geometry = q_cmesh->Reference();
    int dim = geometry->Dimension();
    
    if (!geometry) {
        DebugStop();
    }
   // UniformRefinement(geometry, 1);
    TPZCompMesh *s_cmesh = new TPZCompMesh(geometry);
    
    geometry->ResetReference();
    q_cmesh->LoadReferences();
//    int nel = q_cmesh->NElements();
//    for (int iel =0; iel<nel; iel++) {
//        TPZCompEl *cel = q_cmesh->Element(iel);
//        if (!cel) {
//            continue;
//        }
//        TPZGeoEl *gel = cel->Reference();
//        int nsides = gel->NSides();
//        TPZGeoElSide gelside(gel,nsides-1);
//        TPZGeoElSide neigh = gelside.Neighbour();
//        while (neigh!=gelside) {
//            int mat = neigh.Element()->MaterialId();
//            if (mat == 11 && neigh.Element()->Dimension()==0) {
//                neigh.Element()->SetMaterialId(8);
//            }
//            neigh = neigh.Neighbour();
//        }
//
//    }
    

    std::ofstream filecheck("flux_cmesh_for_transport.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geometry, filecheck);
    
    std::set<int> allmat_ids = {1,2,6,7,8};
    
    for(auto frac_set: fractures)
    {
        for(auto matid:frac_set.m_id)
        {
            allmat_ids.insert(matid);
        }
    }
    
    std::set<int> bcmat_ids = {3,4,5,130,140,150,230,240,250,-1942};
    
    int nstate = 1;
    TPZVec<STATE> sol(1,0.0);
    
    /// Inserting the materials
    for (auto mat_id: allmat_ids) {
        TPZL2Projection * volume = new TPZL2Projection(mat_id,0,nstate, sol);
        s_cmesh->InsertMaterialObject(volume);
    }
    
    TPZMaterial * material = s_cmesh->FindMaterial(1);
    if (!material) {
        DebugStop();
    }
    TPZFMatrix<STATE> val1(1,1,0);
    TPZFMatrix<STATE> val2(1,1,1);
    
    int typ = 0; // inlet
    /// Inserting the materials
    for (auto mat_id: bcmat_ids) {
        TPZMaterial * bc = material->CreateBC(material, mat_id, typ, val1, val2);
        s_cmesh->InsertMaterialObject(bc);
    }
    geometry->ResetReference();
    
    s_cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    int s_order = 0;
    s_cmesh->SetDefaultOrder(s_order);
    for (auto cel : q_cmesh->ElementVec()) {
        
        if (!cel) {
            DebugStop();
        }
        
        TPZGeoEl * gel = cel->Reference();
        int gel_index = gel->Index();
        int mat_id = gel->MaterialId();
        
        if (allmat_ids.find(mat_id) != allmat_ids.end()) {
            TPZGeoEl * gel = geometry->Element(gel_index);
            CreateTransportElement(s_order,s_cmesh, gel, false);
        }
        if (bcmat_ids.find(mat_id) != bcmat_ids.end()) {
            TPZGeoEl * gel = geometry->Element(gel_index);
            CreateTransportElement(s_order,s_cmesh, gel, true);
        }
    }
    
    
    
    s_cmesh->SetDimModel(0);
    s_cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    int ind_0 = -1;
    for (int i = 0; i<fractures.size(); i++) {
        if (fractures[i].m_dim == 0) {
            ind_0 = i;
            break;
        }
    }
    if(ind_0 >= 0)
    {
        /// Create point element
        for (auto cel : cmesh->ElementVec()) {
            if (!cel) {
                continue;
            }
            
            TPZGeoEl * gel = cel->Reference();
            if (!gel) {
                DebugStop();
            }
            
            if (gel->Dimension() != 0) {
                continue;
            }
            int matid = gel->MaterialId();
            if(bcmat_ids.find(matid) != bcmat_ids.end())
            {
                continue;
            }
            if(allmat_ids.find(matid) == allmat_ids.end())
            {
                continue;
            }
            
            TPZMultiphysicsInterfaceElement * int_face = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
            if (int_face) {
                continue;
            }
            
            int n_connect = cel->NConnects();
            bool matid_found = fractures[ind_0].m_id.find(gel->MaterialId())  != fractures[ind_0].m_id.end();
            if (n_connect !=1 && matid_found) {
                DebugStop();
            }
            else if(n_connect == 1 && !matid_found)
            {
                std::cout << "I dont understand matid should be included in the fracture data structure\n";
            }
            else if (n_connect == 0)
            {
                DebugStop();
            }
            TPZConnect & c = cel->Connect(0);
            
            if (c.NElConnected() > 2 && matid_found) {
                CreateTransportElement(s_order,s_cmesh, gel, false);
            }
            
            
        }
    }
    geometry->ResetReference();
    s_cmesh->SetDimModel(dim);
    s_cmesh->InitializeBlock();
    
    return s_cmesh;
}

void CreateTransportElement(int p_order, TPZCompMesh *cmesh, TPZGeoEl *gel, bool is_BC){
    int64_t cel_index;
    int dimension = gel->Dimension();
    cmesh->SetDimModel(dimension);
    if (is_BC) {
        cmesh->SetDimModel(dimension+1);
    }
    TPZCompEl * cel = cmesh->ApproxSpace().CreateCompEl(gel, *cmesh, cel_index);
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
    TPZCompElDisc *intelDisc = dynamic_cast<TPZCompElDisc *> (cel);
    if (intel){
        intel->PRefine(p_order);
    } else if (intelDisc) {
        intelDisc->SetDegree(p_order);
        intelDisc->SetTrueUseQsiEta();
    } else {
        DebugStop();
    }
    gel->ResetReference();
}
TPZAnalysis * CreateTransportAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh, true);
    
    if (sim_data.UsePardisoQ) {
        
        TPZSpStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        analysis->SetSolver(step);
        
        return analysis;
    }else{
        
        TPZSkylineNSymStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        analysis->SetSolver(step);
        analysis->SetStructuralMatrix(matrix);
        return analysis;
        
    }
    
}
TPZFMatrix<STATE> TimeForward(TPZAnalysis * tracer_analysis, int & n_steps, REAL & dt, TPZFMatrix<STATE> & M_diag){
    
    TPZMultiphysicsCompMesh * cmesh_transport = dynamic_cast<TPZMultiphysicsCompMesh *>(tracer_analysis->Mesh());
    
    if (!cmesh_transport) {
        DebugStop();
    }
    
    TPZManVector<TPZCompMesh *,3> meshtrvec = cmesh_transport->MeshVector();
    
    /// Compute mass matrix M.
    TPZAutoPointer<TPZMatrix<STATE> > M;
    TPZFMatrix<REAL> F_inlet;
    {
        
        bool mass_matrix_Q = true;
        std::set<int> volumetric_mat_ids = {1,2,6,7,8};
        
        for (auto mat_id: volumetric_mat_ids) {
            TPZMaterial * mat = cmesh_transport->FindMaterial(mat_id);
            TPZTracerFlow * volume = dynamic_cast<TPZTracerFlow * >(mat);
            if (!volume) {
                continue;
            }
            volume->SetMassMatrixAssembly(mass_matrix_Q);
        }
        
        std::cout << "Computing Mass Matrix." << std::endl;
        tracer_analysis->Assemble();
        std::cout << "Mass Matrix is computed." << std::endl;
        M = tracer_analysis->Solver().Matrix()->Clone();
    }
    int n_rows = M->Rows();
    M_diag.Resize(n_rows,1);
    for (int64_t i = 0; i < n_rows; i++) {
        M_diag(i,0) = M->Get(i, i);
    }
    int64_t n_eq = tracer_analysis->Mesh()->NEquations();
    TPZFMatrix<STATE> saturations(n_eq,n_steps);
    
    {
        bool mass_matrix_Q = false;
        std::set<int> volumetric_mat_ids = {1,2,6,7,8};
        
        for (auto mat_id: volumetric_mat_ids) {
            TPZMaterial * mat = cmesh_transport->FindMaterial(mat_id);
            TPZTracerFlow * volume = dynamic_cast<TPZTracerFlow * >(mat);
            if (!volume) {
                continue;
            }
            volume->SetTimeStep(dt);
            volume->SetMassMatrixAssembly(mass_matrix_Q);
        }
        
        std::cout << "Computing transport operator K = M + T, and F_inlet " << std::endl;
        tracer_analysis->Assemble();
        F_inlet = tracer_analysis->Rhs();
    }
    
    /// Time evolution
    std::string file_reservoir("transport.vtk");
    
    {
        int div = 0;
        TPZFMatrix<REAL> s_n(n_eq,1,0.0);
        TPZFMatrix<REAL> last_state_mass(n_eq,1,0.0);
        TPZFMatrix<REAL> s_np1;
        
        for (int it = 0; it < n_steps; it++) {
            
            for (int64_t i = 0; i < n_eq; i++) {
                last_state_mass(i,0) = M_diag(i,0)*s_n(i,0);
            }
            
            tracer_analysis->Rhs() = F_inlet - last_state_mass;
            tracer_analysis->Rhs() *= -1.0;
            
            tracer_analysis->Solve(); /// (LU decomposition)
            s_np1 = tracer_analysis->Solution();
            tracer_analysis->LoadSolution(s_np1);
            
            /// postprocess ...
            TPZStack<std::string,10> scalnames, vecnames;
            scalnames.Push("Sw");
            scalnames.Push("So");
            
            std::map<int,int> volumetric_ids;
            volumetric_ids.insert(std::make_pair(1, 3));
            
            std::map<int,int> fracture_ids;
            fracture_ids.insert(std::make_pair(6, 2));
            
            std::map<int,int> fracture_intersections_ids;
            fracture_intersections_ids.insert(std::make_pair(7, 1));
            
            for (auto data: volumetric_ids) {
                TPZMaterial * mat = cmesh_transport->FindMaterial(data.first);
                TPZTracerFlow * volume = dynamic_cast<TPZTracerFlow * >(mat);
                if (!volume) {
                    DebugStop();
                }
                volume->SetDimension(data.second);
            }
            
                int div = 0;
                std::set<int> mat_id_3D;
                mat_id_3D.insert(1);
                mat_id_3D.insert(2);
                std::string file_reservoir("cube_s.vtk");
                tracer_analysis->DefineGraphMesh(2,mat_id_3D,scalnames,vecnames,file_reservoir);
                tracer_analysis->PostProcess(div,2);
                
//                std::set<int> mat_id_2D;
//                mat_id_2D.insert(6);
//                std::string file_frac("fracture_s.vtk");
//                tracer_analysis->DefineGraphMesh(1,mat_id_2D,scalnames,vecnames,file_frac);
//                tracer_analysis->PostProcess(div,1);
            
            
            
            // configuring next time step
            s_n = s_np1;
            for (int64_t i = 0; i < n_eq; i++) {
                saturations(i,it) = cmesh_transport->Solution()(i,0);
            }
            
        }
        
        
    }
    return saturations;
}
void UniformRefinement(TPZGeoMesh * geometry, int h_level)
 {
    
    TPZManVector<TPZGeoEl*> sons;
    for(int i=0; i < h_level; i++)
    {
        int64_t nels = geometry->NElements();
        for(int64_t elem = 0; elem < nels; elem++)
        {
            TPZGeoEl * gel = geometry->ElementVec()[elem];
            if(!gel || gel->HasSubElement())
                continue;
            gel->Divide(sons);
        }
    }
    geometry->ResetConnectivities();
    geometry->BuildConnectivity();
}
void Case_1(){
    
    
    /// Execution logger.
    std::ofstream log_file("case_1_results.txt");
    
    int h_level = 0;
    
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.IsHybrid=true;
    sim.n_threads = 24;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0e-5);
    sim.porosities.push_back(0.25);
    
    /// not used but inserted
    sim.omega_ids.push_back(2);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0e-6);
    sim.porosities.push_back(0.2);
    
    /// C inlet value
    sim.c_inlet = 0.01;
    
    int bc_inlet  = 3;
    int bc_outlet = 4;
    int bc_non_flux = 5;
    
    sim.gamma_ids.push_back(bc_inlet);
    sim.gamma_dim.push_back(3);
    sim.gamma_ids.push_back(bc_outlet);
    sim.gamma_dim.push_back(3);
    sim.gamma_ids.push_back(bc_non_flux);
    sim.gamma_dim.push_back(3);
    
    int bc_type_D = 0;    //    D = 0;
    int bc_type_N = 1;    //    N = 1;
    REAL p_inlet  = 4.0;
    REAL p_outlet = 1.0;
    REAL qn       = 0.0;
    
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_N);
    
    sim.vals.push_back(p_inlet);
    sim.vals.push_back(p_outlet);
    sim.vals.push_back(qn);
    
    /// Defining DFN boundary data (id,bc_type,data)
    std::vector<std::tuple<int,int,REAL>> bc_ids_2d;
    bc_ids_2d.push_back(std::make_tuple(bc_inlet,bc_type_D,p_inlet));
    bc_ids_2d.push_back(std::make_tuple(bc_outlet,bc_type_D,p_outlet));
    bc_ids_2d.push_back(std::make_tuple(bc_non_flux,bc_type_N,qn));
    
    //    int bc_1d_inlet  = 130;
    //    int bc_1d_outlet = 140;
    int bc_1d_non_flux = 150;
    
    //    int bc_0d_inlet  = 230;
    //    int bc_0d_outlet = 240;
    //    int bc_0d_non_flux = 250;
    
    std::map<int,int> bc_ids_1d_map;
    //    bc_ids_1d_map.insert(std::make_pair(bc_inlet,bc_1d_inlet));
    //    bc_ids_1d_map.insert(std::make_pair(bc_outlet,bc_1d_outlet));
    bc_ids_1d_map.insert(std::make_pair(bc_non_flux,bc_1d_non_flux));
    
    //    std::map<int,int> bc_ids_0d_map;
    //    bc_ids_0d_map.insert(std::make_pair(bc_inlet,bc_0d_inlet));
    //    bc_ids_0d_map.insert(std::make_pair(bc_outlet,bc_0d_outlet));
    //    bc_ids_0d_map.insert(std::make_pair(bc_non_flux,bc_0d_non_flux));
    
    /// Defining DFN data
    TPZStack<TFracture> fracture_data;
    TFracture fracture;
    fracture.m_id.insert(6);
    fracture.m_dim              = 2;
    fracture.m_kappa_normal     = 20.0;
    fracture.m_kappa_tangential = 1.0e-3;
    fracture.m_d_opening        = 1.0e-2;
    fracture.m_porosity         = 0.4;
    fracture_data.push_back(fracture);
    fracture.m_id.clear();
    
    /// Benchmarks Material ID convention
    /// 1 and 2 for 3D matrix
    /// 3,4, and 5 for 2D matrix boundaries (3 -> inlet, 4 -> outlet, 5 -> impervious)
    /// 6 fractures
    /// 7 fractures intersections
    /// 8 crossing intersections
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag(4); // From 0D to 3D
    dim_name_and_physical_tag[3]["RockMatrix_1"] = 1;
    dim_name_and_physical_tag[3]["RockMatrix_2"] = 2;
    dim_name_and_physical_tag[2]["BCInlet"] = 3;
    dim_name_and_physical_tag[2]["BCOutlet"] = 4;
    dim_name_and_physical_tag[2]["BCImpervious"] = 5;
    dim_name_and_physical_tag[2]["Fractures"] = 6;
    dim_name_and_physical_tag[1]["FracturesIntersections"] = 7;
    dim_name_and_physical_tag[0]["CrossingIntresections"] = 8;
    
    TPZGmshReader Geometry;
    std::string source_dir = SOURCE_DIR;
        std::string file_gmsh = source_dir + "/meshes/Case_1/case_1.msh";
    //    std::string file_gmsh = source_dir + "/meshes/Case_1/case_1_1k.msh";
    //    std::string file_gmsh = source_dir + "/meshes/Case_1/case_1_10k.msh";
    //    std::string file_gmsh = source_dir + "/meshes/Case_1/case_1_100k.msh";
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    std::string version("4.1");
    Geometry.SetFormatVersion(version);
    Geometry.SetDimNamePhysical(dim_name_and_physical_tag);
    gmesh = Geometry.GeometricGmshMesh(file_gmsh.c_str());
    Geometry.PrintPartitionSummary(std::cout);
    
    UniformRefinement(gmesh, h_level);
    
#ifdef PZDEBUG
    std::ofstream file("geometry_case_1_base.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    std::ofstream file_txt("geometry_case_1_base.txt");
    gmesh->Print(file_txt);
#endif
    
    
    int p_order = 2;
    TPZVec<TPZCompMesh *> meshvec;
    //TPZCompMesh *cmixedmesh = NULL;
  //  cmixedmesh = MPCMeshMixed(gmesh, p_order, sim, meshvec);
    //#ifdef PZDEBUG
    //    std::ofstream filemixed("mixed_cmesh.txt");
    //    cmixedmesh->Print(filemixed);
    //#endif
    
    CreateDFN dfn;
    int flux_trace_id=10;
    int lagrange_id=11;
    int flux_resistivity_id=13;
    int mp_nterface_id =12;
    
    dfn.SetFractureData(fracture_data);
    dfn.SetReservoirBoundaryData(bc_ids_2d);
    dfn.SetMapReservoirBCToDFNBC1DIds(bc_ids_1d_map);
 //   dfn.SetMapReservoirBCToDFNBC0DIds(bc_ids_0d_map);
    
    //  UniformRefinement(gmesh, 1);
    TPZMultiphysicsCompMesh *mp_initial_cmesh = MPCMeshMixed(gmesh, p_order, sim, meshvec);
    dfn.SetPeriferalMaterialIds(flux_trace_id, lagrange_id, mp_nterface_id,flux_resistivity_id);
    TPZCompMesh * dfn_hybrid_cmesh = dfn.CreateDFNCmesh(mp_initial_cmesh);
    
    TPZMultiphysicsCompMesh * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(dfn_hybrid_cmesh);
    
    TPZCompMesh *s_cmesh = CreateTransportMesh(mp_cmesh, fracture_data);
    std::ofstream printsat("s_cmesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(s_cmesh, printsat);
    
    std::ofstream printsat2("s_cmesh.txt");
    s_cmesh->Print(printsat2);
    
    TPZManVector<TPZCompMesh *,3> meshtrvec(3);
    meshtrvec[0] = meshvec[0];
    meshtrvec[1] = meshvec[1];
    meshtrvec[2] = s_cmesh;
    
    TPZMultiphysicsCompMesh *cmesh_transport = MPTransportMesh(mp_cmesh, fracture_data, sim, meshtrvec);
    TPZAnalysis * tracer_analysis = CreateTransportAnalysis(cmesh_transport, sim);

    bool solve_dfn_problem_Q = true;
    if (solve_dfn_problem_Q) {
        TPZManVector<TPZCompMesh * > mesh_vec = mp_cmesh->MeshVector();
        
        std::ofstream printpress("pres.txt");
        mesh_vec[1]->Print(printpress);
        //     TPZVTKGeoMesh::PrintCMeshVTK(mesh_vec[1], printpress);
        
        std::ofstream printflux("flux.txt");
        mesh_vec[0]->Print(printflux);
        //      TPZVTKGeoMesh::PrintCMeshVTK(mesh_vec[1], printflux);
        
        std::cout << "Condensing DFN equations." << std::endl;
        std::cout << "DFN neq before condensation = " << mp_cmesh->NEquations() << std::endl;
        dfn.GroupElements(mp_cmesh);
        std::cout << "DFN neq = " << mp_cmesh->NEquations() << std::endl;
        
#ifdef PZDEBUG
        {
            
            std::ofstream file_hybrid_condensed("Hybrid_mixed_condensed.txt");
            mp_cmesh->ComputeNodElCon();
            mp_cmesh->Print(file_hybrid_condensed);
        }
#endif
        
        std::cout << "DFN equations are condensed." << std::endl;
        
        TPZAnalysis *an = CreateAnalysis(mp_cmesh, sim);
        std::cout << "Assembly DFN problem neq = " << mp_cmesh->NEquations() << std::endl;
        an->Assemble();
        std::cout << "Assembly for DFN complete." << std::endl;
        
        std::cout << "Solving DFN problem." << std::endl;
        an->Solve();
        std::cout << "DFN problem solved." << std::endl;
        mp_cmesh->LoadSolutionFromMultiPhysics();
        
        
        TPZStack<std::string,10> scalnames, vecnames;
        vecnames.Push("q");
        scalnames.Push("p");
        
        int div = 0;
        std::set<int> mat_id_3D;
        mat_id_3D.insert(1);
        mat_id_3D.insert(2);
        std::string file_reservoir("cube.vtk");
        an->DefineGraphMesh(3,mat_id_3D,scalnames,vecnames,file_reservoir);
        an->PostProcess(div,3);
        
        std::set<int> mat_id_2D;
        mat_id_2D.insert(6);
        std::string file_frac("fracture.vtk");
        an->DefineGraphMesh(2,mat_id_2D,scalnames,vecnames,file_frac);
        an->PostProcess(div,2);
        
        std::set<int> mat_id_1D;
        mat_id_1D.insert(7);
        std::string file_frac_intersections("fracture_intersections.vtk");
        an->DefineGraphMesh(1,mat_id_1D,scalnames,vecnames,file_frac_intersections);
      //  an->PostProcess(div,1);
    }
  
    
    int n_steps = 100;
    REAL dt     = 1.0e7;
    TPZFMatrix<STATE> M_diag;
    TPZFMatrix<STATE> saturations = TimeForward(tracer_analysis, n_steps, dt, M_diag);
    
    
    return;
}

void TesisDamner(){
    
    //Lectura de la malha
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag(4); // From 0D to 3D
    dim_name_and_physical_tag[2]["RockMatrix_1"] = 1;
    dim_name_and_physical_tag[2]["RockMatrix_2"] = 2;
    dim_name_and_physical_tag[1]["BCInlet"] = 3;
    dim_name_and_physical_tag[1]["BCOutlet"] = 4;
    dim_name_and_physical_tag[1]["BCImpervious"] = 5;
    dim_name_and_physical_tag[1]["Fractures"] = 6;
    dim_name_and_physical_tag[1]["Fractures_2"] = 6;
    dim_name_and_physical_tag[1]["Fractures_3"] = 6;
    dim_name_and_physical_tag[0]["FracturesIntersections"] = 7;
    dim_name_and_physical_tag[0]["CrossingIntresections"] = 8;
    
    TPZGmshReader Geometry;
    std::string source_dir = SOURCE_DIR;
    std::string file_gmsh = source_dir + "/meshes/Damner/test_damner2.msh";
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    std::string version("4.1");
    Geometry.SetFormatVersion(version);
    
    Geometry.SetDimNamePhysical(dim_name_and_physical_tag);
    gmesh = Geometry.GeometricGmshMesh(file_gmsh.c_str());
    Geometry.PrintPartitionSummary(std::cout);
    
   
    
#ifdef PZDEBUG
    std::ofstream file("geometry_cube_base.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    std::ofstream file_txt("geometry_cube_base.txt");
    gmesh->Print(file_txt);
    gmesh->SetDimension(2);
#endif
  //  UniformRefinement(gmesh, 1);
    //Configuracion de Materiales
    
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.IsHybrid=true;
    sim.n_threads = 1;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(2);
    sim.permeabilities.push_back(1.0);
    sim.porosities.push_back(1.0);
    
    
    /// C inlet value
    sim.c_inlet = 1.0;
    
    int bc_inlet  = 3;
    int bc_outlet = 4;
    int bc_non_flux = 5;
    
    sim.gamma_ids.push_back(bc_inlet);
    sim.gamma_dim.push_back(2);
    sim.gamma_ids.push_back(bc_outlet);
    sim.gamma_dim.push_back(2);
    sim.gamma_ids.push_back(bc_non_flux);
    sim.gamma_dim.push_back(2);
    
    int bc_type_D = 0;    //    D = 0;
    int bc_type_N = 1;    //    N = 1;
    REAL p_inlet  = 2.0;
    REAL p_outlet = 1.0;
    REAL qn       = 0.0;
    
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_N);
    
    sim.vals.push_back(p_inlet);
    sim.vals.push_back(p_outlet);
    sim.vals.push_back(qn);
    TPZVec<TPZCompMesh *> meshvec;
    
    /// Defining DFN data
    TPZStack<TFracture> fracture_data;
    TFracture fracture;
    fracture.m_id.insert(6);
    fracture.m_dim              = 1;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 0.5;
    fracture_data.push_back(fracture);
    fracture.m_id.clear();
    fracture.m_id.insert(7);
    fracture.m_dim              = 1;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 1.0;
    fracture_data.push_back(fracture);
    fracture.m_id.clear();
    fracture.m_id.insert(8);
    fracture.m_dim              = 0;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 1.0;
    fracture_data.push_back(fracture);
    
    /// Defining DFN boundary data (id,bc_type,data)
    std::vector<std::tuple<int,int,REAL>> bc_ids_2d;
    bc_ids_2d.push_back(std::make_tuple(bc_inlet,bc_type_D,p_inlet));
    bc_ids_2d.push_back(std::make_tuple(bc_outlet,bc_type_D,p_outlet));
    bc_ids_2d.push_back(std::make_tuple(bc_non_flux,bc_type_N,qn));
    
    int bc_1d_inlet  = 130;
    int bc_1d_outlet = 140;
    int bc_1d_non_flux = 150;
    
    int bc_0d_inlet  = 230;
    int bc_0d_outlet = 240;
    int bc_0d_non_flux = 250;
    
    std::map<int,int> bc_ids_1d_map;
    bc_ids_1d_map.insert(std::make_pair(bc_inlet,bc_1d_inlet));
    bc_ids_1d_map.insert(std::make_pair(bc_outlet,bc_1d_outlet));
    bc_ids_1d_map.insert(std::make_pair(bc_non_flux,bc_1d_non_flux));
    
    std::map<int,int> bc_ids_0d_map;
    bc_ids_0d_map.insert(std::make_pair(bc_inlet,bc_0d_inlet));
    bc_ids_0d_map.insert(std::make_pair(bc_outlet,bc_0d_outlet));
    bc_ids_0d_map.insert(std::make_pair(bc_non_flux,bc_0d_non_flux));
    //
    CreateDFN dfn;
    int flux_trace_id=10;
    int lagrange_id=11;
    int flux_resistivity_id=13;
    int mp_nterface_id =12;
    
    dfn.SetFractureData(fracture_data);
    dfn.SetReservoirBoundaryData(bc_ids_2d);
    dfn.SetMapReservoirBCToDFNBC1DIds(bc_ids_1d_map);
    dfn.SetMapReservoirBCToDFNBC0DIds(bc_ids_0d_map);
    
  //  UniformRefinement(gmesh, 1);
    TPZMultiphysicsCompMesh *mp_initial_cmesh = MPCMeshMixed(gmesh, 1, sim, meshvec);
    dfn.SetPeriferalMaterialIds(flux_trace_id, lagrange_id, mp_nterface_id,flux_resistivity_id);
    TPZCompMesh * dfn_hybrid_cmesh = dfn.CreateDFNCmesh(mp_initial_cmesh);
    
    TPZMultiphysicsCompMesh * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(dfn_hybrid_cmesh);
    
    
    TPZCompMesh *s_cmesh = CreateTransportMesh(mp_cmesh, fracture_data);
    std::ofstream printsat("s_cmesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(s_cmesh, printsat);

    std::ofstream printsat2("s_cmesh.txt");
    s_cmesh->Print(printsat2);
    
    TPZManVector<TPZCompMesh *,3> meshtrvec(3);
    meshtrvec[0] = meshvec[0];
    meshtrvec[1] = meshvec[1];
    meshtrvec[2] = s_cmesh;

    TPZMultiphysicsCompMesh *cmesh_transport = MPTransportMesh(mp_cmesh, fracture_data, sim, meshtrvec);
    TPZAnalysis * tracer_analysis = CreateTransportAnalysis(cmesh_transport, sim);
    
    
    bool solve_dfn_problem_Q = true;
    if (solve_dfn_problem_Q) {
        TPZManVector<TPZCompMesh * > mesh_vec = mp_cmesh->MeshVector();
        
        std::ofstream printpress("pres.vtk");
        //      mesh_vec[1]->Print(printpress);
        TPZVTKGeoMesh::PrintCMeshVTK(mesh_vec[1], printpress);
        
        std::ofstream printflux("flux.vtk");
        //      mesh_vec[0]->Print(printflux);
        TPZVTKGeoMesh::PrintCMeshVTK(mesh_vec[0], printflux);
        
        std::cout << "Condensing DFN equations." << std::endl;
        std::cout << "DFN neq before condensation = " << mp_cmesh->NEquations() << std::endl;
        //  dfn.GroupElements(mp_cmesh);
        std::cout << "DFN neq = " << mp_cmesh->NEquations() << std::endl;
        
#ifdef PZDEBUG
        {
            
            //            std::ofstream file_hybrid_condensed("Hybrid_mixed_condensed.txt");
            //            mp_cmesh->ComputeNodElCon();
            //            mp_cmesh->Print(file_hybrid_condensed);
        }
#endif
        
        std::cout << "DFN equations are condensed." << std::endl;
        
        TPZAnalysis *an = CreateAnalysis(mp_cmesh, sim);
        std::cout << "Assembly DFN problem neq = " << mp_cmesh->NEquations() << std::endl;
        an->Assemble();
        std::cout << "Assembly for DFN complete." << std::endl;
        
        std::cout << "Solving DFN problem." << std::endl;
        an->Solve();
        std::cout << "DFN problem solved." << std::endl;
        mp_cmesh->LoadSolutionFromMultiPhysics();
        
        
        TPZStack<std::string,10> scalnames, vecnames;
        vecnames.Push("q");
        scalnames.Push("p");
        
        int div = 0;
        std::set<int> mat_id_3D;
        mat_id_3D.insert(1);
        //        mat_id_3D.insert(2);
        std::string file_reservoir("caso1_damner.vtk");
        an->DefineGraphMesh(2,mat_id_3D,scalnames,vecnames,file_reservoir);
        an->PostProcess(div,2);
        
//        std::set<int> mat_id_2D;
//        mat_id_2D.insert(6);
//        std::string file_frac("fracture.vtk");
//        an->DefineGraphMesh(2,mat_id_2D,scalnames,vecnames,file_frac);
//        an->PostProcess(div,2);
//
//        std::set<int> mat_id_1D;
//        mat_id_1D.insert(6);
//        std::string file_frac_intersections("fracture_intersections.vtk");
//        an->DefineGraphMesh(1,mat_id_1D,scalnames,vecnames,file_frac_intersections);
//        an->PostProcess(div,1);
    }
    
    int n_steps = 10;
    REAL dt     = 10.0;
    TPZFMatrix<STATE> M_diag;
    TPZFMatrix<STATE> saturations = TimeForward(tracer_analysis, n_steps, dt, M_diag);

    
    
    
    return;
}
