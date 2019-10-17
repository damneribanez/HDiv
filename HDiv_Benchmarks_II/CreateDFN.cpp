    
#include "CreateDFN.h"


#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "pzvec.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr dfn_logger(Logger::getLogger("DFN"));
#endif


/// Default constructor
CreateDFN::CreateDFN() : TPZHybridizeHDiv(){
    
    m_fracture_ids.clear();
    m_fracture_data.resize(0);
    m_geometry = NULL;
}

/// Default destructor
CreateDFN::~CreateDFN(){
    
}

void CreateDFN::CreateInterfaceElements(int target_dim, int interface_id, TPZCompMesh *cmesh, TPZVec<TPZCompMesh *> & mesh_vec) {
    
    TPZCompMesh * pressure_cmesh = mesh_vec[1];
    TPZVec<TPZCompEl *> celpressure(pressure_cmesh->NElements(), 0);
    for (int64_t el = 0; el < pressure_cmesh->NElements(); el++) {
        TPZCompEl *cel = pressure_cmesh->Element(el);
        if (cel && cel->Reference() && cel->Reference()->Dimension() == target_dim - 1) {
            celpressure[el] = cel;
        }
    }
    cmesh->Reference()->ResetReference();
    for (auto &cel :cmesh->ElementVec()) {
        if(cel->Reference()->Dimension() == target_dim -1){
            cel->LoadElementReference();
        }
    }
    
    TPZManVector<int64_t,2> left_mesh_indexes(1,1), right_mesh_indexes(1,0);
    
    for (auto cel : celpressure) {
        if (!cel) continue;
        TPZStack<TPZCompElSide> celstack;
        TPZGeoEl *gel = cel->Reference();
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        
        TPZCompElSide celside = gelside.Reference();
        
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        int count = 0;
        for (auto &celstackside : celstack) {
            if (celstackside.Reference().Element()->Dimension() == target_dim - 1) {
                TPZGeoElBC gbc(gelside, interface_id);
                
                TPZCompEl * right_cel = celstackside.Element();
                int n_connects = right_cel->NConnects();
                if (n_connects == 1) {
                    int64_t index;
                    
                    TPZMultiphysicsInterfaceElement * mp_interface_el = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside, celstackside);
                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                    count++;
                }
                
            }
        }
        if (count == 1)
        {
            TPZCompElSide clarge = gelside.LowerLevelCompElementList2(false);
            if(!clarge) DebugStop();
            TPZGeoElSide glarge = clarge.Reference();
            if (glarge.Element()->Dimension() == target_dim) {
                TPZGeoElSide neighbour = glarge.Neighbour();
                while(neighbour != glarge)
                {
                    if (neighbour.Element()->Dimension() < target_dim) {
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
                if(neighbour == glarge) DebugStop();
                glarge = neighbour;
            }
            clarge = glarge.Reference();
            if(!clarge) DebugStop();
            TPZGeoElBC gbc(gelside, interface_id);
            
            int64_t index;
            TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside, clarge);
            mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
            count++;
        }
    }
    pressure_cmesh->InitializeBlock();
}


void CreateDFN::SetFractureData(TPZStack<TFracture> & fracture_data){
    m_fracture_data = fracture_data;
    /// Fill m_fracture_ids
    for (auto fracture: m_fracture_data) {
        m_fracture_ids[fracture.m_dim] = fracture.m_id;
    }
}




/// insert the given material ids in the fluxmesh.
// insert the impervious material to seal the fractures
// put the lagrange matid in the pressure mesh
void CreateDFN::InsertMaterials(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & flux_resistivity_id){
    
    int n_state = 1;
    
    TPZMultiphysicsCompMesh  * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(cmesh);
    if (!mp_cmesh) {
        DebugStop();
    }
    TPZManVector<TPZCompMesh *, 3> dfn_mixed_mesh_vec = mp_cmesh->MeshVector();
    
    /// Insert LM and flux_trace materials
    {
        TPZVec<STATE> sol(1,0);
        TPZCompMesh * flux_cmesh = dfn_mixed_mesh_vec[0];
        if (!flux_cmesh) {
            DebugStop();
        }
        if (!flux_cmesh->FindMaterial(flux_trace_id)) {
            auto flux_trace_mat = new TPZL2Projection(flux_trace_id, target_dim, n_state, sol);
            flux_trace_mat->SetScaleFactor(0.0);
            flux_cmesh->InsertMaterialObject(flux_trace_mat);
        }
        
        if (!flux_cmesh->FindMaterial(flux_resistivity_id)) {
            auto flux_trace_mat = new TPZL2Projection(flux_resistivity_id, target_dim, n_state, sol);
            flux_trace_mat->SetScaleFactor(0.0);
            flux_cmesh->InsertMaterialObject(flux_trace_mat);
        }
        
        
        TPZCompMesh * pressure_cmesh = dfn_mixed_mesh_vec[1];
        if (!pressure_cmesh) {
            DebugStop();
        }
        if (!pressure_cmesh->FindMaterial(lagrange_id)) {
            auto lagrange_mult_mat = new TPZL2Projection(lagrange_id, target_dim, n_state, sol);
            lagrange_mult_mat->SetScaleFactor(0.0);
            pressure_cmesh->InsertMaterialObject(lagrange_mult_mat);
        }
    }
    
}



TPZVec<int> & CreateDFN::ExtractActiveApproxSpaces(TPZCompMesh * cmesh){
    
    TPZMultiphysicsCompMesh  * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(cmesh);
    if (!mp_cmesh) {
        DebugStop();
    }
    
    TPZVec<int> & active_approx_spaces = mp_cmesh->GetActiveApproximationSpaces();
    return active_approx_spaces;
}

TPZVec<TPZCompMesh *> & CreateDFN::GetMeshVector(TPZCompMesh * cmesh){
    TPZMultiphysicsCompMesh  * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(cmesh);
    if (!mp_cmesh) {
        DebugStop();
    }
    return mp_cmesh->MeshVector();
}

TPZCompMesh * CreateDFN::DuplicateMultiphysicsCMeshMaterials(TPZCompMesh * cmesh){
    
    TPZMultiphysicsCompMesh  * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(cmesh);
    if (!mp_cmesh) {
        DebugStop();
    }
    
    TPZGeoMesh *gmesh = cmesh->Reference();
    TPZMultiphysicsCompMesh * dfn_hybrid_cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->CopyMaterials(*dfn_hybrid_cmesh);
    return dfn_hybrid_cmesh;
}

void CreateDFN::CleanUpMultiphysicsCMesh(TPZCompMesh * cmesh){
    cmesh->CleanUp();
}

void CreateDFN::InsertMaterialsForHibridization(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & flux_resistivity_id, int & mp_nterface_id){
    
    int n_state = 1;
    TPZVec<STATE> sol(1,0);
    if (!cmesh) {
        DebugStop();
    }
    if (!cmesh->FindMaterial(flux_trace_id)) {
        
        auto flux_trace_mat = new TPZL2Projection(flux_trace_id, target_dim, n_state, sol);
        flux_trace_mat->SetScaleFactor(0.0);
        cmesh->InsertMaterialObject(flux_trace_mat);
        
    }
    
    /// TODO:: Use this structure for get normal permeability
    REAL kappa_normal = m_fracture_data[0].m_kappa_normal;
    if (!cmesh->FindMaterial(flux_resistivity_id)) {
        
        auto normal_flux_mat = new TPZNormalDarcyFlow(flux_resistivity_id, target_dim);
        normal_flux_mat->SetKappaNormal(kappa_normal);
        cmesh->InsertMaterialObject(normal_flux_mat);
        
    }
    
    if (!cmesh->FindMaterial(lagrange_id)) {
        
        auto lagrange_mult_mat = new TPZL2Projection(lagrange_id, target_dim, n_state, sol);
        lagrange_mult_mat->SetScaleFactor(0.0);
        cmesh->InsertMaterialObject(lagrange_mult_mat);
        
    }
    
    if (!cmesh->FindMaterial(mp_nterface_id)) {
        auto lagrange_multiplier = new TPZLagrangeMultiplier(mp_nterface_id, target_dim - 1, n_state);
        lagrange_multiplier->SetMultiplier(1.0);
        cmesh->InsertMaterialObject(lagrange_multiplier);
    }
    
}

void CreateDFN::BuildFracturesAproxSpace(int p_order, int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & flux_resistivity_id, int & mp_nterface_id) {
    
    TPZMultiphysicsCompMesh  * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(cmesh);
    if (!mp_cmesh) {
        DebugStop();
    }
    TPZManVector<TPZCompMesh *, 3> dfn_mixed_mesh_vec = mp_cmesh->MeshVector();
    
    
    /// Change lagrange multipliers identifiers to fractures identifiers
    {
        
        /// insert material objects in pressure space
        {
            TPZCompMesh * pressure_cmesh = dfn_mixed_mesh_vec[1];
            pressure_cmesh->Reference()->ResetReference();
            pressure_cmesh->LoadReferences();
            int n_state = 1;
            TPZVec<STATE> sol(1,0);
            for(auto fracture : m_fracture_data){
                std::set<int> fracture_ids = fracture.m_id;
                for (auto matid: fracture_ids)
                {
                    auto fracture_material = new TPZL2Projection(matid, target_dim, n_state, sol);
                    fracture_material->SetScaleFactor(0);
                    if (!pressure_cmesh->FindMaterial(matid)) {
                        pressure_cmesh->InsertMaterialObject(fracture_material);
                    }
                    else
                    {
                        delete fracture_material;
                    }
                }
            }
            
        }
        
        TPZCompMesh * pressure_cmesh = dfn_mixed_mesh_vec[1];
        pressure_cmesh->Reference()->ResetReference();
        pressure_cmesh->LoadReferences();
        TPZGeoMesh  * gmesh = pressure_cmesh->Reference();
      
        // switch the material id of pressure elements if they are neighbour of fracture elements
        for (auto gel : gmesh->ElementVec()) {
            
            if (!gel) {
                DebugStop();
            }
            
            if (gel->MaterialId() != lagrange_id) {
                continue;
            }
            // cel is a pressure element
            TPZCompEl * cel = gel->Reference();
            if (!cel) {
                continue;
            }
            
            // if the material id is a lagrange material id, look if there is no geometric element
            // neighbour of the pressure element that has a material id of a fracture element
            if (gel->MaterialId() == lagrange_id) {
                int side = gel->NSides() - 1;
                int geldim = gel->Dimension();
                TPZGeoElSide gel_side(gel,side);
                TPZStack<TPZGeoElSide> gelsides;
                gel_side.AllNeighbours(gelsides);
                
                bool fracture_detected_Q;
                for (auto iside : gelsides) {
                    TPZGeoEl * neigh_gel = iside.Element();
                    int gel_mat_id = neigh_gel->MaterialId();
                    std::set<int>::iterator it = m_fracture_ids[geldim].find(gel_mat_id);
                    fracture_detected_Q = it != m_fracture_ids[geldim].end();
                    // We found a fracture element and are switching the geometric reference of the computational pressure element
                    // There should be a flux element associated with this geometric element as well!
                    if (fracture_detected_Q && neigh_gel->Dimension() == geldim) {
                        cel->SetReference(neigh_gel->Index());
                    }
                    
                }
                
                /// change the material id of HDivBound to resistivity id
                if (fracture_detected_Q) { /// Change Hdivbound to Normal flux operator
                    for (auto iside : gelsides) {
                        TPZGeoEl * neigh_gel = iside.Element();
                        int gel_mat_id = neigh_gel->MaterialId();
                        if (gel_mat_id != flux_trace_id) {
                            continue;
                        }
                        neigh_gel->SetMaterialId(flux_resistivity_id);
                    }
                }
                
            }
            
        }
        
    }
    
    
    /// create boundary elements for the 2d and 1d fractures
    
    TPZCompMesh * flux_cmesh = dfn_mixed_mesh_vec[0];
    flux_cmesh->Reference()->ResetReference();
    int n_state = 1;
    TPZVec<STATE> sol(1,0);
    // insert the material objects if needed
    std::set<int> fracture_flux_material_set;
    int c = 0;
    for(auto fracture : m_fracture_data){
        if (fracture.m_dim != target_dim) {
            continue;
        }
        std::set<int> fracture_id = fracture.m_id;
        for (auto matid : fracture_id)
        {
            auto fracture_material = new TPZL2Projection(matid, target_dim, n_state, sol);
            fracture_material->SetScaleFactor(0.0);
            if (!flux_cmesh->FindMaterial(matid)) {
                flux_cmesh->InsertMaterialObject(fracture_material);
                fracture_flux_material_set.insert(matid);
            }
            else
            {
                delete fracture_material;
            }
            auto frac_mat = flux_cmesh->FindMaterial(matid);
            if (c == 0) {
                
                TPZFMatrix<STATE> val1(1,1,0.0), val2(1,1,0.0);
                for (auto chunk : m_bc_ids_2d) {
                    
                    int bc_mat_id(std::get<0>(chunk));
                    int bc_type(std::get<1>(chunk));
                    REAL bc_value(std::get<2>(chunk));
                    
                    int bc_fracture_id;
                    if (target_dim == 2) {
                        std::cout<<"entré"<<std::endl;
                        bc_fracture_id = m_bc_ids_1d[bc_mat_id];
                        val2(0,0) = bc_value;
                        if (!flux_cmesh->FindMaterial(bc_fracture_id)) {
                            auto bc = frac_mat->CreateBC(frac_mat, bc_fracture_id, bc_type, val1, val2);
                            flux_cmesh->InsertMaterialObject(bc);
                            fracture_flux_material_set.insert(bc_fracture_id);
                        }
                    }else if (target_dim == 1){
                        bc_fracture_id = m_bc_ids_0d[bc_mat_id];
                        val2(0,0) = bc_value;
                        if (!flux_cmesh->FindMaterial(bc_fracture_id)) {
                            auto bc = frac_mat->CreateBC(frac_mat, bc_fracture_id, bc_type, val1, val2);
                            flux_cmesh->InsertMaterialObject(bc);
                            fracture_flux_material_set.insert(bc_fracture_id);
                        }
                    }
                    
                }
            }
            c = 1;
        }
    }
    
    
    std::set<int> bc_indexes_2d, bc_indexes_1d, bc_indexes_0d;
    for (auto chunk : m_bc_ids_2d) {
        int mat_id(std::get<0>(chunk));
        bc_indexes_2d.insert(mat_id);
    }
    
    /// Insert fractures intersections with 3D bc
    TPZGeoMesh * geometry = flux_cmesh->Reference();
    geometry->SetDimension(target_dim);
    switch (target_dim) {
        case 2:
        {
            for (auto chunk : m_bc_ids_1d) {
                int mat_id = chunk.second;
                bc_indexes_1d.insert(mat_id);
            }
            CreateFractureBCGeoElements(target_dim, geometry, bc_indexes_2d, bc_indexes_1d);
        }
            break;
        case 1:
        {
            for (auto chunk : m_bc_ids_0d) {
                int mat_id = chunk.second;
                bc_indexes_0d.insert(mat_id);
            }
            CreateFractureBCGeoElements(target_dim, geometry, bc_indexes_2d, bc_indexes_0d);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
#ifdef PZDEBUG
    std::ofstream geo_file("geometry_with_bcs.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geometry, geo_file, true);
#endif
    
    /// create the target_dim fracture elements (including the boundary conditions)
    LoadReferencesByDimension(flux_cmesh, target_dim);
    flux_cmesh->SetDimModel(target_dim);
    flux_cmesh->SetDefaultOrder(p_order);
    flux_cmesh->SetAllCreateFunctionsHDiv();
    flux_cmesh->AutoBuild(fracture_flux_material_set);

    /// this step makes the flux mesh "not computable". If all sideorients are one, the mesh without hybridizing is inconsistent
    int sideorient = 1;
    for (auto gel : flux_cmesh->Reference()->ElementVec()) {
        if (!gel->Reference()) {
            continue;
        }
        
        int n_sides = gel->NSides();
        
        for (int is = 0 ; is < n_sides; is++) {
            if (gel->SideDimension(is) != target_dim-1) {
                continue;
            }
            
            TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement *>(gel->Reference());
            
            if (!intel) {
                DebugStop();
            }
            
            intel->SetSideOrient(is, sideorient);
        }
    }
    
    flux_cmesh->InitializeBlock();
    flux_cmesh->SetDimModel(target_dim); ///  for coherence
    
}

/// create boundary elements of dimension target_dim-1 which lay on the boundary of the three dimensional mesh
void CreateDFN::CreateFractureBCGeoElements(int target_dim, TPZGeoMesh * gmesh, std::set<int> bc_indexes, std::set<int> bc_frac_indexes){
    
    
    int dim = gmesh->Dimension();
    std::map<std::pair<int,int>,std::pair<int,int>> surf_to_surf_side_indexes;
    
    for (auto gel : gmesh->ElementVec()) {
        
        if (!gel) continue;
        // TODO dim will always be 3 : so we only consider 2d elements????
        if (gel->Dimension() != dim - 1) continue;
        if (gel->HasSubElement()) {
            continue;
        }
        
        bool quad_gel_Q = gel->Type() == EQuadrilateral;
        bool trin_gel_Q = gel->Type() == ETriangle;
        bool line_gel_Q = gel->Type() == EOned;
        
        std::vector<int> sides;
        if (quad_gel_Q) {
            if (target_dim == 2) {
                sides = {4,5,6,7};
            }else{
                sides = {0,1,2,3};
            }
            
        }else if (trin_gel_Q){
            
            if (target_dim == 2) {
                sides = {3,4,5};
            }else{
                sides = {0,1,2};
            }
        }else if (line_gel_Q){
            sides = {0,1};
        }
        
        int bc_mat_id = gel->MaterialId();
        if(bc_indexes.find(bc_mat_id) ==  bc_indexes.end())
        {
            continue;
        }
        
        // here we have a geometric element that corresponds to a 3d boundary (2 dimensional)
        for (auto side: sides) {
            TPZStack<TPZGeoElSide> all_neigh;
            TPZGeoElSide gelside(gel, side);
            // for each side of the boundary element, look for all neighbours
            gelside.AllNeighbours(all_neigh);
            std::set<int> surfaces;
            surfaces.insert(gel->Index());
            
            // foundface indicates whether there is another geometric element with the same material id
            bool foundface = false;
            // foundbc indicates we found an element with material id of fracture bcs
            bool foundbc = false;
            // bcmatid will be equal the boundary material of a neighbouring face
            int bcmatid = 0;
            // foundfrac indicates there is a neighbouring geometric element with a fracture material id
            bool foundfrac = false;
            for (auto gel_side : all_neigh) {
                bool is_bc_member_Q = gel_side.Element()->MaterialId() == bc_mat_id;
                if (is_bc_member_Q) {
                    foundface = true;
                    // this statement is redundant???
                    bcmatid = gel_side.Element()->MaterialId();
                }
                if(bc_frac_indexes.find(gel_side.Element()->MaterialId()) != bc_frac_indexes.end())
                {
                    foundbc = true;
                    
                }
                
                if(m_fracture_ids[target_dim].find(gel_side.Element()->MaterialId()) != m_fracture_ids[target_dim].end())
                {
                    if(gel_side.Element()->Dimension() == target_dim) {
                        foundfrac = true;
                    }
                }
            }
            bool candidate_to_create_bc_Q = foundface == false || foundbc == true || foundfrac == false;
            if(candidate_to_create_bc_Q) continue;
            for (auto gel_side : all_neigh) {
                if(m_fracture_ids[target_dim].find(gel_side.Element()->MaterialId()) != m_fracture_ids[target_dim].end())
                {
                    if (target_dim == 2) {
                        int bc_frac_mat_id = m_bc_ids_1d[bcmatid];
                        TPZGeoElBC gbc(gelside,bc_frac_mat_id);
                    }else if (target_dim == 1){
                        int bc_frac_mat_id = m_bc_ids_0d[bcmatid];
                        TPZGeoElBC gbc(gelside,bc_frac_mat_id);
                    }
                    
                    break;
                }
            }
        }
    }
    
}

void CreateDFN::InsertMaterialsForMixedOperatorOnFractures(int target_dim, TPZCompMesh * cmesh){
    
    int c = 0;
    for(auto fracture : m_fracture_data){
        for(auto matid:fracture.m_id)
        {
            if (!cmesh->FindMaterial(matid))
            {
                if (fracture.m_dim != target_dim)
                {
                    continue;
                }
                if(fracture.m_dim != 0)
                {
                    auto fracture_material = new TPZMixedDarcyFlow(matid, target_dim);
                    fracture_material->SetPermeability(fracture.m_kappa_tangential);
                    cmesh->InsertMaterialObject(fracture_material);
                    if (c == 0)
                    {
                        
                        TPZFMatrix<STATE> val1(1,1,0.0), val2(1,1,0.0);
                        for (auto chunk : m_bc_ids_2d) {
                            
                            int bc_mat_id(std::get<0>(chunk));
                            int bc_type(std::get<1>(chunk));
                            REAL bc_value(std::get<2>(chunk));
                            
                            int bc_fracture_id;
                            if (target_dim == 2) {
                                bc_fracture_id = m_bc_ids_1d[bc_mat_id];
                                val2(0,0) = bc_value;
                                if (!cmesh->FindMaterial(bc_fracture_id)) {
                                    auto bc = fracture_material->CreateBC(fracture_material, bc_fracture_id, bc_type, val1, val2);
                                    cmesh->InsertMaterialObject(bc);
                                }
                            }else if (target_dim == 1){
                                bc_fracture_id = m_bc_ids_0d[bc_mat_id];
                                val2(0,0) = bc_value;
                                if (!cmesh->FindMaterial(bc_fracture_id)) {
                                    auto bc = fracture_material->CreateBC(fracture_material, bc_fracture_id, bc_type, val1, val2);
                                    cmesh->InsertMaterialObject(bc);
                                }
                            }
                            
                        }
                        c++;
                    }
                }
                else // the fracture dimension is zero, create an L2 material
                {
                    TPZVec<STATE> sol(1,0.);
                    int n_state = 1;
                    auto lagrange_mat = new TPZL2Projection(matid, target_dim, n_state, sol);
                    lagrange_mat->SetScaleFactor(0.0);
                    cmesh->InsertMaterialObject(lagrange_mat);
                    
                }
            }
        }
    }
    
}

void CreateDFN::LoadReferencesByDimension(TPZCompMesh * flux_cmesh, int dim){
    
#ifdef PZDEBUG
    if (!flux_cmesh) {
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geometry = flux_cmesh->Reference();
    geometry->ResetReference();
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    for (auto cel: flux_cmesh->ElementVec()) {
        
        if (!cel) {
            continue;
        }
        
        TPZGeoEl * gel = cel->Reference();
        if (!gel) {
            continue;
        }
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        bool check_vols = gel->Dimension() == dim && cel->NConnects() > 1; ///  Volumetric elements
        bool check_bc = gel->Dimension() == dim - 1 && cel->NConnects() == 1; ///  Boundary elements
        
        if (check_vols || check_bc) {
            cel->LoadElementReference();
            
#ifdef LOG4CXX
            if(dfn_logger->isDebugEnabled())
            {
                std::stringstream sout;
                sout << " LoadReferencesByDimension:: ";
                sout << " gel index =  " << gel->Index();
                sout << " matid = " << gel->MaterialId();
                sout << " connect indexes =  ";
                int n_connects = cel->NConnects();
                for (int i =0; i < n_connects; i++) {
                    sout << cel->ConnectIndex(i) << " ";
                }
                LOGPZ_DEBUG(dfn_logger,sout.str())
            }
#endif
            
        }
    }
}



void CreateDFN::BuildMultiphysicsCMesh(int dim, TPZCompMesh * hybrid_cmesh, TPZVec<int> & approx_spaces, TPZManVector<TPZCompMesh *, 3> mesh_vec){
    
    TPZMultiphysicsCompMesh  * mp_hybrid_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(hybrid_cmesh);
    if (!mp_hybrid_cmesh) {
        DebugStop();
    }
    
    mp_hybrid_cmesh->SetDimModel(dim);
    mp_hybrid_cmesh->BuildMultiphysicsSpace(approx_spaces, mesh_vec);
}

/// group and condense the elements
void CreateDFN::GroupElements(TPZMultiphysicsCompMesh *cmesh)
{
    TPZVec<int64_t> ConnectGroup(cmesh->NConnects(),-1);
    TPZVec<int64_t> ElementGroup(cmesh->NElements());
    int64_t nelem = cmesh->NElements();
    int64_t igr = 0;
    // each element that has a flux associated is the nucleus of a group
    for (int64_t el=0; el<nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZMultiphysicsElement *mpel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if(!mpel) continue;
        TPZCompEl *atom_cel_flux = mpel->Element(0);
        if(!atom_cel_flux)
        {
            continue;
        }
        int ncon_flux = atom_cel_flux->NConnects();
        if(ncon_flux == 1)
        {
            continue;
        }
        for(int ic=0; ic<ncon_flux; ic++)
        {
            int64_t conindex = cel->ConnectIndex(ic);
            ConnectGroup[conindex] = igr;
        }
        ElementGroup[el] = igr;
        igr++;
    }
    int64_t numgroups = igr;
    // an element that shares a connect with the nucleus element will belong to the group
    for (int64_t el=0; el<nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        int ncon = cel->NConnects();
        int64_t elgr = -1;
        for(int ic=0; ic<ncon; ic++)
        {
            int64_t conindex = cel->ConnectIndex(ic);
            if(ConnectGroup[conindex] != -1) elgr = ConnectGroup[conindex];
        }
        ElementGroup[el] = elgr;
    }
    TPZVec<TPZElementGroup *> groups(numgroups,0);
    for (igr = 0; igr < numgroups; igr++) {
        int64_t index;
        TPZElementGroup *group = new TPZElementGroup(*cmesh,index);
        groups[igr] = group;
    }
    for (int64_t el=0; el<nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(ElementGroup[el] != -1)
        {
            int64_t group_index = ElementGroup[el];
            TPZElementGroup *group = groups[group_index];
            group->AddElement(cel);
        }
    }
    cmesh->ComputeNodElCon();
    for (int igr = 0; igr < numgroups; igr++) {
        TPZElementGroup *group = groups[igr];
        bool keep_matrix = false;
        new TPZCondensedCompEl(group, keep_matrix);
    }
    cmesh->CleanUpUnconnectedNodes();
}

TPZCompMesh * CreateDFN::CreateDFNCmesh(TPZMultiphysicsCompMesh * mp_initial_cmesh){
    TPZManVector<TPZCompMesh *> mp_mesh_vector = mp_initial_cmesh->MeshVector();
    TPZCompMesh * q_mesh = mp_mesh_vector[0];
    TPZCompMesh * p_mesh = mp_mesh_vector[1];
    q_mesh->SetDimModel(2);
    q_mesh->Reference()->SetDimension(2);
    int order = 1;
    std::ofstream file1("p_sides.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(p_mesh, file1);
    std::ofstream file2("q_sides.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(q_mesh, file2);
    TPZHybridizeHDiv hybrid;
    hybrid.SetPeriferalMaterialIds(fHDivWrapMatid, fLagrangeInterface, fInterfaceMatid);
    hybrid.ComputeNState(mp_initial_cmesh->MeshVector());
    hybrid.InsertPeriferalMaterialObjects(mp_initial_cmesh->MeshVector());
    hybrid.HybridizeInternalSides(mp_initial_cmesh->MeshVector());
    std::ofstream file3("p_Hybrid_internal_sides.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(p_mesh, file3);
    std::ofstream file4("q_Hybrid_internal_sides.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(q_mesh, file4);
 

    for (int dim = 2; dim > 1; dim--) {
        std::cout<<"dim: "<<dim<<std::endl;
         hybrid.HybridizeInternalSides(mp_initial_cmesh->MeshVector());
        BuildFracturesAproxSpace(order, dim, mp_initial_cmesh, fHDivWrapMatid, fLagrangeInterface, fflux_resistivity_id, fInterfaceMatid);
        hybrid.HybridizeInternalSides(mp_initial_cmesh->MeshVector());

        std::ofstream file5("geometry_new.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(p_mesh->Reference(), file5);
        std::ofstream file3("p_Hybrid_int_build.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(p_mesh, file3);
        std::ofstream file4("q_Hybrid_int_build.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(q_mesh, file4);
        std::ofstream file6("pressure.txt");
        p_mesh->Print(file6);

    }
    /// Computational multiphysics mesh reconstruction
    TPZVec<int> & active_approx_spaces = ExtractActiveApproxSpaces(mp_initial_cmesh);
    TPZCompMesh * dfn_hybrid_cmesh = DuplicateMultiphysicsCMeshMaterials(mp_initial_cmesh);
    CleanUpMultiphysicsCMesh(mp_initial_cmesh);
    
    int target_dim=2;
    int fractures_dim =1;
    int fractures_intersections_dim =0;
    InsertMaterialsForHibridization(target_dim, dfn_hybrid_cmesh, fHDivWrapMatid, fLagrangeInterface, fflux_resistivity_id, fInterfaceMatid);
    InsertMaterialsForMixedOperatorOnFractures(fractures_dim,dfn_hybrid_cmesh);
    InsertMaterialsForMixedOperatorOnFractures(fractures_intersections_dim,dfn_hybrid_cmesh);
    InsertMaterialsForMixedOperatorOnFractures(0,dfn_hybrid_cmesh);
    
    
    TPZManVector<TPZMaterial *> zeroDimMat(m_fracture_ids[0].size());
    int count = 0;
    for (auto matid:m_fracture_ids[0]) {
        zeroDimMat[count++] = dfn_hybrid_cmesh->FindMaterial(matid);
        dfn_hybrid_cmesh->MaterialVec().erase(matid);
    }
    BuildMultiphysicsCMesh(target_dim,dfn_hybrid_cmesh,active_approx_spaces,mp_mesh_vector);
   
    
    for (auto matptr:zeroDimMat) {
        dfn_hybrid_cmesh->InsertMaterialObject(matptr);
    }
  
    CreateAllInterfaceElements(fInterfaceMatid, dfn_hybrid_cmesh, mp_mesh_vector);
    dfn_hybrid_cmesh->ComputeNodElCon();
    CheckODElements(dfn_hybrid_cmesh);

    dfn_hybrid_cmesh->InitializeBlock();
    return dfn_hybrid_cmesh;
}
/// make a hybrid mesh from a H(div) multiphysics mesh
void *CreateDFN::Hybridize(TPZMultiphysicsCompMesh *multiphysics, bool group_elements, double Lagrange_term_multiplier)
{
    
    TPZHybridizeHDiv hybrid;
    hybrid.SetPeriferalMaterialIds(fHDivWrapMatid, fLagrangeInterface, fInterfaceMatid);
    hybrid.ComputeNState(multiphysics->MeshVector());
    hybrid.InsertPeriferalMaterialObjects(multiphysics->MeshVector());
    hybrid.HybridizeInternalSides(multiphysics->MeshVector());
    
}
void CreateDFN::CheckODElements(TPZCompMesh *dfn_hybrid_cmesh){

        int64_t nel = dfn_hybrid_cmesh->NElements();
        for(int64_t el=0; el<nel; el++)
        {
            TPZCompEl *cel = dfn_hybrid_cmesh->Element(el);
            TPZMultiphysicsElement *mpcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if(!mpcel) continue;
            if(!mpcel->Element(0) && mpcel->Element(1))
            {
                
                TPZConnect &c = mpcel->Connect(0);
                int nelcon = c.NElConnected();
                if(nelcon <= 3)
                {
                    continue;
                }
                std::cout << "this element has a pressure but no flux\n";
                
                TPZGeoEl *gel = cel->Reference();
                // find a neighbour with matid = 8
                TPZGeoEl *gel8 = 0;
                int matid_target = 0;
                TPZGeoElSide gelside(gel,gel->NSides()-1);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside)
                {
                    if(m_fracture_ids[0].find(neighbour.Element()->MaterialId()) != m_fracture_ids[0].end())
                    {
                        gel8 = neighbour.Element();
                        matid_target = neighbour.Element()->MaterialId();
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
                /// if there is no element 8 found, issue a warning because the mesher should have created one
                if(!gel8)
                {
                    std::cout << "no element found at the intersection of fractures\n";
                    int matid = *m_fracture_ids[0].begin();
                    std::cout << "changin matid to " << matid << std::endl;
                    gel->SetMaterialId(matid);
                }
                else
                {
                    gel8->SetMaterialId(0);
                    gel8->RemoveConnectivities();
                    gel->SetMaterialId(matid_target);
                }
            }
        }
}
void CreateDFN::CreateAllInterfaceElements(int interface_id, TPZCompMesh *cmesh, TPZVec<TPZCompMesh *> & mesh_vec){
        CreateInterfaceElements(3, fInterfaceMatid, cmesh, mesh_vec);
        CreateInterfaceElements(2, fInterfaceMatid, cmesh, mesh_vec);
        CreateInterfaceElements(1, fInterfaceMatid, cmesh, mesh_vec);
}
