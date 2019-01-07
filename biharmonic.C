// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include files needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "mesh_tools.h"
#include "mesh_subdiv_support.h"
#include "mesh_modification.h"
#include "mesh_refinement.h"
#include "vtk_io.h"
#include "equation_systems.h"
#include "linear_implicit_system.h"

// Define the Finite Element object.
#include "fe.h"
#include "dof_map.h"
#include "dof_map_subdiv.h"
#include "utility.h"
#include "getpot.h"

// Define Gauss quadrature rules.
#include "quadrature_gauss.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "petsc_vector.h"
#include "petsc_matrix.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "elem.h"

#include "petsc_macro.h"
EXTERN_C_FOR_PETSC_BEGIN
# include <petscts.h>
EXTERN_C_FOR_PETSC_END


typedef struct {
    RealVectorValue a[2];  //Tangent vectors (2D surface)
    DenseMatrix<Real> gcov;  //Covariant metric tensor components
    DenseMatrix<Real> gcon;  //Contravariant metric tensor components
    RealVectorValue acovd[2][2];  //Covariant derivatives of tangent vectors a_(a,b)
    PetscReal g_abc[2][2][2];  //Covariant derivatives of metric tensor g_(ab,c)
    std::vector<Real> LaplaceSphi;
    std::vector<RealVectorValue> GradSphi;
    RealVectorValue n;  //Surface normal (not normalized)
    PetscScalar J1;    //Length of surface normal = Jacobian of mapping from parametric to physical space
    PetscScalar invgdet;  // 1/det(gcov)   => needed for construction of contravariant metric
    std::vector<DenseMatrix<Real> > dphicov2;
} SurfaceData;

class CachedElement
{
public:
    CachedElement() {initialized = false; element_fe=NULL; }
    ~CachedElement() { if (initialized) { delete element_fe; element_fe = NULL; delete element_qrule; element_qrule = NULL; initialized=false;  }}
	const Elem* elem;
	FEBase* element_fe;
    QBase *element_qrule;
	std::vector<unsigned int> dof_indices_u;
	bool initialized;
	bool is_ghost;
    std::vector<SurfaceData> SD;
};

std::vector<CachedElement> ELEMENT_CACHE;

double CHR(int m, int i, int j,  SurfaceData &SD)
{
    int k; double s = 0.0;
    for (k=0; k<2; k++)
    {
        s += 0.5*SD.gcon(m,k)*(SD.g_abc[k][i][j] + SD.g_abc[k][j][i] - SD.g_abc[i][j][k]);
    }
    return s;
}

//CHECKED AND CORRECT
void setup_surface_data(int qp, SurfaceData& SD, FEBase* element_fe)
{
	START_LOG("setup_surface_metrics", "ShellElement");
    
    SD.gcov.resize(2,2);
    SD.gcon.resize(2,2);
    int n_dofs = element_fe->get_phi().size();
    SD.LaplaceSphi.resize(n_dofs);
    SD.GradSphi.resize(n_dofs);
    SD.dphicov2.resize(n_dofs);

	const std::vector<std::vector<RealGradient> >& dphi = element_fe->get_dphi();
	const std::vector<std::vector<RealTensor> >& d2phi = element_fe->get_d2phi();
	const std::vector< RealGradient >& dxyzdxi = element_fe->get_dxyzdxi();
	const std::vector< RealGradient >& dxyzdeta = element_fe->get_dxyzdeta();
	const std::vector< RealGradient >& d2xyzdxi2 = element_fe->get_d2xyzdxi2();
	const std::vector< RealGradient >& d2xyzdeta2 = element_fe->get_d2xyzdeta2();
	const std::vector< RealGradient >& d2xyzdxideta = element_fe->get_d2xyzdxideta();
	RealGradient ddispdxi, ddispdeta, d2dispdxi2, d2dispdeta2, d2dispdxideta;
	Real det;
		
	SD.a[0] = dxyzdxi[qp];
	SD.a[1] = dxyzdeta[qp];
	SD.n = SD.a[0].cross(SD.a[1]);
	//SD.invgdet = sqrt(ar.acov[2]*ar.acov[2]);
	//SD.invgdet = 1/(SD.n.size());
	
	SD.acovd[0][0] = d2xyzdxi2[qp];
	SD.acovd[0][1] = SD.acovd[1][0] = d2xyzdxideta[qp];
	SD.acovd[1][1] = d2xyzdeta2[qp];
	
	// Covariant metric tensor:
    SD.gcov(0,0) = SD.a[0]*SD.a[0];
    SD.gcov(1,0) = SD.gcov(0,1) = SD.a[0]*SD.a[1];
    SD.gcov(1,1) = SD.a[1]*SD.a[1];
	//Contravariant metric tensor:
	det = (SD.gcov(0,0)*SD.gcov(1,1) - SD.gcov(1,0)*SD.gcov(1,0));
	SD.invgdet = 1./det;
    SD.gcon(0,0) = SD.invgdet*SD.gcov(1,1);
    SD.gcon(1,0) = SD.gcon(0,1) = -SD.invgdet*SD.gcov(0,1);
    SD.gcon(1,1) = SD.invgdet*SD.gcov(0,0);
	
    //covariant derivatives of metric tensor g_(ab,c)
    SD.g_abc[0][0][0] = 2. * SD.a[0]*SD.acovd[0][0];
    SD.g_abc[0][0][1] = 2. * SD.a[0]*SD.acovd[0][1];
    SD.g_abc[0][1][0] = SD.g_abc[1][0][0] = SD.a[1]*SD.acovd[0][0] + SD.acovd[1][0]*SD.a[0];
    SD.g_abc[0][1][1] = SD.g_abc[1][0][1] = SD.a[1]*SD.acovd[0][1] + SD.acovd[1][1]*SD.a[0];
    SD.g_abc[1][1][0] = 2. * SD.a[1]*SD.acovd[1][0];
    SD.g_abc[1][1][1] = 2. * SD.a[1]*SD.acovd[1][1];
    
    for (int a=0; a<n_dofs; a++) {
        double laplace = 0.0;
        //Contract with contravariant metric:
        int j,k;
	DenseMatrix<Real> cov2;
	cov2.resize(2,2);
        for (j=0; j<2; j++)
        {
            for (k=0; k<2; k++) {
               // laplace += SD.gcon(j,k)*N2[a][j][k] - SD.gcon(j,k)*(CHR(0,j,k,SD.g_abc,SD.gcon)*N1[a][0]
               //                                                       + CHR(1,j,k,SD.g_abc,SD.gcon)*N1[a][1]);
                laplace += SD.gcon(j,k)*d2phi[a][qp](j,k) -SD.gcon(j,k)*(CHR(0,j,k,SD)*dphi[a][qp](0)
                                                                        + CHR(1,j,k,SD)*dphi[a][qp](1));
        	cov2(j,k) = d2phi[a][qp](j,k) - (CHR(0,j,k,SD)*dphi[a][qp](0) + CHR(1,j,k,SD)*dphi[a][qp](1)); 
	    }
        }
        SD.dphicov2[a] = cov2;  //Second covariant derivatives of basis function a
        SD.LaplaceSphi[a]=laplace;
        SD.GradSphi[a] = SD.gcon(0,0)*dphi[a][qp](0)*SD.a[0]+SD.gcon(0,1)*dphi[a][qp](1)*SD.a[0]
        + SD.gcon(1,0)*dphi[a][qp](0)*SD.a[1]+SD.gcon(1,1)*dphi[a][qp](1)*SD.a[1];
        //  PetscPrintf(PETSC_COMM_SELF,"LaplaceSN=%g\n",laplace);
    }

}

void print_surface_data(SurfaceData& SD)
{
    std::cout<<"Surface vectors:"<<std::endl;
    std::cout<<"("<<SD.a[0](0)<<", "<<SD.a[0](1)<<", "<<SD.a[0](2)<<"), ("<<SD.a[1](0)<<", "<<SD.a[1](1)<<", "<<SD.a[1](2)<<")\n";
    std::cout<<"Surface con. metric:"<<std::endl;
    std::cout<<SD.gcon(0,0)<<"\t"<<SD.gcon(0,1)<<"\n"<<SD.gcon(0,1)<<"\t"<<SD.gcon(1,1)<<"\n";
    std::cout<<"Surface cov. metric:"<<std::endl;
    std::cout<<SD.gcov(0,0)<<"\t"<<SD.gcov(0,1)<<"\n"<<SD.gcov(0,1)<<"\t"<<SD.gcov(1,1)<<"\n";
    std::cout<<"Surface cov. metric:"<<std::endl;
    std::cout<<SD.gcov(0,0)<<"\t"<<SD.gcov(0,1)<<"\n"<<SD.gcov(0,1)<<"\t"<<SD.gcov(1,1)<<"\n";
    
}

void init_cached_elements(std::vector<CachedElement>& cache, LinearImplicitSystem &sys)
{
    EquationSystems& es = sys.get_equation_systems();
    const MeshBase& mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    const unsigned int u_var = sys.variable_number ("u");

    FEType fe_type = sys.variable_type(u_var);
    
    cache.resize(mesh.n_elem());
    
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;
        const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
        if (sdelem->is_ghost())
            continue;
        int id=elem->id();
        CachedElement* ce=&(cache[id]);
        ce->elem = elem;
        ce->element_fe = FEBase::build(dim, fe_type).release();
        ce->element_qrule = fe_type.default_quadrature_rule(dim,1).release();
        ce->element_fe->attach_quadrature_rule(ce->element_qrule);
        ce->element_fe->reinit(elem);
        sys.get_dof_map().dof_indices(elem, ce->dof_indices_u, u_var);
       
        int nqp = ce->element_qrule->n_points();
        ce->SD.resize(nqp);
        for (int qp=0; qp<nqp; qp++)
            setup_surface_data(qp, ce->SD[qp], ce->element_fe);
        
        ce->initialized=true;
    }
    
}

void SH_residual(TS ts, double t, Vec u, Vec u_t, Vec u0, double t0, Vec r, void* ctx)
{
    LinearImplicitSystem &sys=
    *(static_cast<LinearImplicitSystem*> (ctx));
    PetscVector<Number> U(u), R(r), U_T(u_t), U0(u0);

    EquationSystems& es = sys.get_equation_systems();
    const MeshBase& mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    const unsigned int u_var = sys.variable_number ("u");
    SurfaceData SD;
   // std::cout<<"Time: "<<t<<", t0: "<<t0<<std::endl;
    R.zero();
    double free_energy=0.;
    //std::cout<<"Current U0:\n";
    //U0.print(std::cout);
    //std::cout<<"Current U:\n";
    //U.print(std::cout);
    //std::cout<<"Current U_T:\n";
    //U_T.print(std::cout);
    
    CachedElement* ce;
    
/*  FEType fe_type = sys.variable_type(u_var);
    AutoPtr<FEBase> fe  (FEBase::build(dim, fe_type));
    AutoPtr<QBase> qrule(fe_type.default_quadrature_rule(dim,1));
    fe->attach_quadrature_rule (qrule.get());
  */
    FEBase* fe;
    
    double par_a = es.parameters.get<Real>("a");
    double par_b = es.parameters.get<Real>("b");
    double par_c = es.parameters.get<Real>("c");
   // double par_f = es.parameters.get<Real>("f");
    
    double gamma0 = es.parameters.get<Real>("gamma0");
    double gamma2 = es.parameters.get<Real>("gamma2");
    double penalty = es.parameters.get<Real>("penalty");
    
    
    DenseVector<Number> Re;
    std::vector<unsigned int> dof_indices_u;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;
        const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
        if (sdelem->is_ghost())
            continue;
        ce = &(ELEMENT_CACHE[elem->id()]);
        dof_indices_u = ce->dof_indices_u;
        const unsigned int n_dofs   = dof_indices_u.size();
        fe = ce->element_fe;
        
        const std::vector<Real>& JxW = fe->get_JxW();
        const std::vector<std::vector<Real> >& phi = fe->get_phi();
        const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
        //const std::vector<std::vector<RealTensor> >& d2phi = fe->get_d2phi();
        
        Re.resize (n_dofs);
        Re.zero();
        //std::cout<<"Number of QP: "<<qrule->n_points()<<std::endl;
        for (unsigned int qp=0; qp<ce->element_qrule->n_points(); qp++)
        {
            SD=ce->SD[qp];
            
            double u0=0;
            double u_t=0;
            RealGradient GradSu;
            GradSu(0)=0; GradSu(1)=0;
            double LaplaceSu=0;
            DenseMatrix<Real> dUcov2;
            dUcov2.resize(2,2);
            double dUdxi[2];
            dUdxi[0]=dUdxi[1]=0.0;
            dUcov2(0,0)=dUcov2(1,1)=dUcov2(0,1)=dUcov2(1,0)=0.0;
            for (unsigned int l=0; l<n_dofs; l++)
            {
                u0 += phi[l][qp]*U(dof_indices_u[l]);
                u_t += phi[l][qp]*U_T(dof_indices_u[l]);
                LaplaceSu += SD.LaplaceSphi[l]*U(dof_indices_u[l]);
                dUdxi[0] += dphi[l][qp](0)*U(dof_indices_u[l]);
                dUdxi[1] += dphi[l][qp](1)*U(dof_indices_u[l]);

                GradSu += SD.GradSphi[l]*U(dof_indices_u[l]);
                dUcov2(0,0) += SD.dphicov2[l](0,0)*U(dof_indices_u[l]);
                dUcov2(0,1) += SD.dphicov2[l](0,1)*U(dof_indices_u[l]);
                dUcov2(1,0) += SD.dphicov2[l](1,0)*U(dof_indices_u[l]);
                dUcov2(1,1) += SD.dphicov2[l](1,1)*U(dof_indices_u[l]);
            }
            double GLTerms = -(par_a*u0 + par_b*u0*u0 + par_c*u0*u0*u0 );
            
            free_energy += JxW[qp]*(0.5*par_a*u0*u0 + 1./3.*par_b*u0*u0*u0 + 0.25*par_c*u0*u0*u0*u0
                                    + gamma0/2.*GradSu*GradSu + gamma2/2.*LaplaceSu*LaplaceSu);

            for (unsigned int a=0; a<n_dofs; a++) {
                //It's + gamma0 due to an additional minus sign from partial integration (and another one for the biharmonic one)
               // RealGradient GradSphi = SD.gcon(0,0)*dphi[a][qp](0)*SD.a[0]+SD.gcon(0,1)*dphi[a][qp](1)*SD.a[0]
               // + SD.gcon(1,0)*dphi[a][qp](0)*SD.a[1]+SD.gcon(1,1)*dphi[a][qp](1)*SD.a[1];
                
                Re(a) += JxW[qp]*(
                               phi[a][qp]*u_t - phi[a][qp]*GLTerms + gamma2*SD.LaplaceSphi[a]*LaplaceSu
                               + gamma0*(GradSu*SD.GradSphi[a])
                               );
            }
        }
       R.add_vector(Re, dof_indices_u);
    }
    el     = mesh.active_local_elements_begin();
    //Dirichlet boundary conditions:
    
    for ( ; el != end_el; ++el)
	{
		Elem* elem = *el;
		const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
        if (sdelem->is_ghost()) {
            for (unsigned int s=0; s<elem->n_sides(); s++)
            {
                if (elem->neighbor(s) == NULL)
                {
                    Node* node = elem->get_node(s);
                    int dofidx = node->dof_number(0,u_var,0);
                    R.add(dofidx, penalty*U(dofidx));
                    int nexts=s+1;
                    if (nexts>2) nexts=0;
                        
                    node = elem->get_node(nexts);
                    dofidx = node->dof_number(0,u_var,0);
                    R.add(dofidx, penalty*U(dofidx));
                }
            }
     
//            Node* node = elem->get_node(0);
//            int dofidx = node->dof_number(0,u_var,0);
//            R.add(dofidx, penalty*U(dofidx));
//            node = elem->get_node(1);
//            dofidx = node->dof_number(0,u_var,0);
//            R.add(dofidx, penalty*U(dofidx));
//            node = elem->get_node(2);
//            dofidx = node->dof_number(0,u_var,0);
//            R.add(dofidx, penalty*U(dofidx));
     
        }
    }
    
    R.close();
//    std::cout<<"u_t:\n";
//    sys.get_vector("u_t").print(std::cout);
//    sys.current_solution.print(std::cout);
   // std::cout<<"\nRHS:\n";
   // R.print(std::cout);
    es.parameters.set<Real>("free_energy") = free_energy;

}

void SH_jacobian(TS ts, double time, Vec u, Vec u_t, double shift, Mat* j, Mat* pc, void* ctx)
{
    LinearImplicitSystem &sys=
    *(static_cast<LinearImplicitSystem*> (ctx));
    PetscVector<Number> U(u), U_T(u_t);
    PetscMatrix<Number> J(*j);
   // PC.attach_dof_map(sys.get_dof_map());
    J.attach_dof_map(sys.get_dof_map());
    
    EquationSystems& es = sys.get_equation_systems();
    const MeshBase& mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    const unsigned int u_var = sys.variable_number ("u");
    SurfaceData SD;
    J.zero();
    
    double gamma0 = es.parameters.get<Real>("gamma0");
    double gamma2 = es.parameters.get<Real>("gamma2");
    double penalty = es.parameters.get<Real>("penalty");
    FEBase* fe;
    CachedElement* ce;
    
    DenseMatrix<Number> Ke;
    std::vector<unsigned int> dof_indices_u;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;
        
        const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
        if (sdelem->is_ghost())
            continue;
        
        ce = &(ELEMENT_CACHE[elem->id()]);
        dof_indices_u = ce->dof_indices_u;
        const unsigned int n_dofs   = dof_indices_u.size();
        fe = ce->element_fe;
        
        const std::vector<Real>& JxW = fe->get_JxW();
        const std::vector<std::vector<Real> >& phi = fe->get_phi();
        const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
        const std::vector<std::vector<RealTensor> >& d2phi = fe->get_d2phi();
        
        Ke.resize (n_dofs,n_dofs);
        
        for (unsigned int qp=0; qp<ce->element_qrule->n_points(); qp++)
        {
            SD=ce->SD[qp];
            for (int a=0; a<n_dofs; a++) {
                
              for (int b=0; b<n_dofs; b++) {
                  //RealGradient GradSphib = SD.gcon(0,0)*dphi[a][qp](0)*SD.a[0]+SD.gcon(0,1)*dphi[a][qp](1)*SD.a[0]
                  //+ SD.gcon(1,0)*dphi[a][qp](0)*SD.a[1]+SD.gcon(1,1)*dphi[a][qp](1)*SD.a[1];
                                  //It's + gamma0 due to an additional minus sign from partial integration (and another one for the biharmonic one)
                  Ke(a,b) += JxW[qp]*(shift*phi[a][qp]*phi[b][qp] + gamma2*(SD.LaplaceSphi[a]*SD.LaplaceSphi[b])
                                       +gamma0*(SD.GradSphi[a]*SD.GradSphi[b])
                                       );
              }
            }
        }
        J.add_matrix(Ke,dof_indices_u);
    }
    J.close();
    el     = mesh.active_local_elements_begin();
    //Dirichlet boundary conditions:
    PetscScalar trace=0.0;
    MatGetTrace(J.mat(), &trace);
    double bcwt=trace/J.n();
    double pn2=penalty*bcwt;
    for (; el != end_el; ++el)
    {
        const Elem* elem = *el;
        const Tri3SD* gh_elem = static_cast<const Tri3SD*> (elem);
        if (!gh_elem->is_ghost())
            continue;
        
        // Find the side which is part of the physical plate boundary,
        // that is, the boundary of the original mesh without ghosts.
        for (unsigned int s=0; s<elem->n_sides(); ++s)
        {
            const Tri3SD* nb_elem = static_cast<const Tri3SD*> (elem->neighbor(s));
            if (nb_elem == NULL || nb_elem->is_ghost())
                continue;
            
            Node* nodes [4]; // n1, n2, n3, n4
            nodes[1] = gh_elem->get_node(s); // n2
            nodes[2] = gh_elem->get_node(MeshTools::Subdiv::next[s]); // n3
            nodes[3] = gh_elem->get_node(MeshTools::Subdiv::prev[s]); // n4
            
            // The node in the interior of the domain, \p n1, is the
            // hardest to find.  Walk along the edges of element \p nb until
            // we have identified it.
            unsigned int n = 0;
            nodes[0] = nb_elem->get_node(0);
            while (nodes[0]->id() == nodes[1]->id() || nodes[0]->id() == nodes[2]->id())
                nodes[0] = nb_elem->get_node(++n);
            
            // The penalty value.  \f$ \frac{1}{\epsilon} \f$
           // const Real penalty = 1.e10;
            
                for (unsigned int n=0; n<4; ++n)
            {
                int u_dof = nodes[n]->dof_number (sys.number(), u_var, 0);
                MatZeroRowsColumns(J.mat(),1,&u_dof,pn2,PETSC_NULL,PETSC_NULL);
                //J.add (u_dof, u_dof, penalty);
            }
        }
    } // end of ghost element loop
    
//    J.close();
  //  J.print(std::cout);
    
}



void random_initial_conditions(EquationSystems& es, const std::string& system_name)
{
    std::cout<<"Setting random initial conditions\n";
    libmesh_assert(system_name == "SH");
    LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("SH");
    MeshBase& mesh = system.get_mesh();
    
    MeshBase::const_node_iterator it = mesh.nodes_begin();
    const MeshBase::const_node_iterator nodes_end = mesh.nodes_end();
    for ( ; it!=nodes_end; ++it)
	{
		Node* node=*it;
		int dofidx =  node->dof_number(0,0,0);
        double posx=(*node)(0);
		double dd=static_cast<Real>(std::rand())/static_cast<Real>(RAND_MAX);
        system.solution->set(dofidx,(dd-0.5)*0.0002);
//		system.solution->set(dofidx,0.001*posx);
	}
    
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    //Dirichlet boundary conditions: Set solution back to zero for ghosts and their neighbouring elements
    for (; el != end_el; ++el)
    {
        const Elem* elem = *el;
        const Tri3SD* gh_elem = static_cast<const Tri3SD*> (elem);
        if (!gh_elem->is_ghost())
            continue;
        
        // Find the side which is part of the physical plate boundary,
        // that is, the boundary of the original mesh without ghosts.
        for (unsigned int s=0; s<elem->n_sides(); ++s)
        {
            const Tri3SD* nb_elem = static_cast<const Tri3SD*> (elem->neighbor(s));
            if (nb_elem == NULL || nb_elem->is_ghost())
                continue;
            
            Node* nodes [4]; // n1, n2, n3, n4
            nodes[1] = gh_elem->get_node(s); // n2
            nodes[2] = gh_elem->get_node(MeshTools::Subdiv::next[s]); // n3
            nodes[3] = gh_elem->get_node(MeshTools::Subdiv::prev[s]); // n4
            
            // The node in the interior of the domain, \p n1, is the
            // hardest to find.  Walk along the edges of element \p nb until
            // we have identified it.
            unsigned int n = 0;
            nodes[0] = nb_elem->get_node(0);
            while (nodes[0]->id() == nodes[1]->id() || nodes[0]->id() == nodes[2]->id())
                nodes[0] = nb_elem->get_node(++n);
            
            for (unsigned int n=0; n<4; ++n)
            {
                const unsigned int u_dof = nodes[n]->dof_number (system.number(), 0, 0);
                system.solution->set(u_dof,0);
            }
        }
    } // end of ghost element loop

}

void planar_stripes_initial_conditions(EquationSystems& es, const std::string& system_name)
{
    std::cout<<"Setting planar stripes initial conditions\n";
    libmesh_assert(system_name == "SH");
    LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("SH");
    MeshBase& mesh = system.get_mesh();
    
    MeshBase::const_node_iterator it = mesh.nodes_begin();
    const MeshBase::const_node_iterator nodes_end = mesh.nodes_end();
    for ( ; it!=nodes_end; ++it)
	{
		Node* node=*it;
		int dofidx =  node->dof_number(0,0,0);
        double posx=(*node)(0);
        double posy=(*node)(1);

        double dd;
     /*   if (fabs(posx)<=1) {
            dd=0.0001*cos(0.5*3.14159*posx);
            system.solution->set(dofidx,dd);
        //		system.solution->set(dofidx,0.001*posx);
        } else if (fabs(posy)<=1) {
            dd=0.0001*cos(0.5*3.14159*posy);
            system.solution->set(dofidx,dd);
        }*/
        if (fabs(posx)<=1 && posy<0) {
            dd=0.0001*cos(0.5*3.14159*posx);
            system.solution->set(dofidx,dd);
            //		system.solution->set(dofidx,0.001*posx);
        }
	}
    
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    //Dirichlet boundary conditions: Set solution back to zero for ghosts and their neighbouring elements
    for (; el != end_el; ++el)
    {
        const Elem* elem = *el;
        const Tri3SD* gh_elem = static_cast<const Tri3SD*> (elem);
        if (!gh_elem->is_ghost())
        continue;
        
        // Find the side which is part of the physical plate boundary,
        // that is, the boundary of the original mesh without ghosts.
        for (unsigned int s=0; s<elem->n_sides(); ++s)
        {
            const Tri3SD* nb_elem = static_cast<const Tri3SD*> (elem->neighbor(s));
            if (nb_elem == NULL || nb_elem->is_ghost())
            continue;
            
            Node* nodes [4]; // n1, n2, n3, n4
            nodes[1] = gh_elem->get_node(s); // n2
            nodes[2] = gh_elem->get_node(MeshTools::Subdiv::next[s]); // n3
            nodes[3] = gh_elem->get_node(MeshTools::Subdiv::prev[s]); // n4
            
            // The node in the interior of the domain, \p n1, is the
            // hardest to find.  Walk along the edges of element \p nb until
            // we have identified it.
            unsigned int n = 0;
            nodes[0] = nb_elem->get_node(0);
            while (nodes[0]->id() == nodes[1]->id() || nodes[0]->id() == nodes[2]->id())
            nodes[0] = nb_elem->get_node(++n);
            
            for (unsigned int n=0; n<4; ++n)
            {
                const unsigned int u_dof = nodes[n]->dof_number (system.number(), 0, 0);
                system.solution->set(u_dof,0);
            }
        }
    } // end of ghost element loop
    
}

void winded_stripes_initial_conditions(EquationSystems& es, const std::string& system_name)
{
    std::cout<<"Setting winded stripes initial conditions\n";
    libmesh_assert(system_name == "SH");
    LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("SH");
    MeshBase& mesh = system.get_mesh();
    double winding_ic_a = es.parameters.get<Real>("winding_ic_a");
    double winding_ic_w = es.parameters.get<Real>("winding_ic_w");
    double winding_ic_my = es.parameters.get<Real>("winding_ic_my");

    MeshBase::const_node_iterator it = mesh.nodes_begin();
    const MeshBase::const_node_iterator nodes_end = mesh.nodes_end();
    for ( ; it!=nodes_end; ++it)
	{
		Node* node=*it;
		int dofidx =  node->dof_number(0,0,0);
        double px=(*node)(0);
        double py=(*node)(1);
        double pz=(*node)(2);
        if (py<winding_ic_my) {
            if (fabs(pz-winding_ic_a*px)/winding_ic_w<1.0)
            {
                double dd=0.0001*cos(0.5*3.14159*(pz-winding_ic_a*px)/winding_ic_w);
                system.solution->set(dofidx,dd);

            }
            
        }
	}
    
}

void project_to_1_sphere(MeshBase& mesh)
{
    MeshBase::const_node_iterator           nd = mesh.nodes_begin();
    const MeshBase::const_node_iterator end_nd = mesh.nodes_end();
    for ( ; nd != end_nd; ++nd)
    {
        Node* node = *nd;
        *node /=  node->size();
    }
    
}

extern "C"
{
    PetscErrorCode sh_petsc_residual (TS ts, PetscReal t, Vec u, Vec u_t, Vec r, void *ctx)
    {
        int ierr=0;
       // std::cout<<"calling residual\n"<<std::endl;
        double tt=static_cast<double>(t);
        LinearImplicitSystem &sys=
        *(static_cast<LinearImplicitSystem*> (ctx));
        Vec u0;
        PetscReal      t0;
        ierr = TSGetTime(ts,&t0);CHKERRQ(ierr);
        ierr = TSGetSolution(ts,&u0);CHKERRQ(ierr);
        SH_residual(ts, tt, u, u_t, u0, static_cast<double>(t0), r, ctx);
        
        
       /*
        PetscVector<Number> U_input(u), R_input(r), U_t_input(u_t);
        std::cout<<"U_t petsc:\n";
         VecView(u_t,PETSC_VIEWER_STDOUT_WORLD);
        std::cout<<"U petsc:\n";
        VecView(u,PETSC_VIEWER_STDOUT_WORLD);
        
        PetscVector<Number>& U_system = *dynamic_cast<PetscVector<Number>*>(sys.solution.get());
        PetscVector<Number>* U0_system = dynamic_cast<PetscVector<Number>*>(&(sys.get_vector("u0")));
//        PetscMatrix<double>* J = dynamic_cast<PetscMatrix<double>*>(&jac);
        Vec u0; //= U0_system->vec();
        PetscVector<Number>& R_system = *dynamic_cast<PetscVector<Number>*>(sys.rhs);
        PetscVector<Number>& U_t_system = *dynamic_cast<PetscVector<Number>*>(&(sys.get_vector("u_t")));
        //Get the last known solution and store it in U0_system:
        ierr = TSGetSolution(ts,&u0);CHKERRQ(ierr);
       // std::cout<<"U0 petsc:\n";
       // VecView(u0,PETSC_VIEWER_STDOUT_WORLD);
        PetscVector<Number> U0(u0);
        U0_system->swap(U0);
        std::cout<<"U0_system:\n";
        U0_system->print(std::cout);
        
        U_input.swap(U_system);
        U_t_input.swap(U_t_system);
        U_input.swap(R_system);
        
        sys.get_dof_map().enforce_constraints_exactly(sys);
        sys.update();
        
        SH_residual(sys,tt);
        
        // Swap back
        R_system.close();
        U_input.swap(U_system);
        U_t_input.swap(U_t_system);
        R_input.swap(R_system);
      //  R_input.close();

      //  std::cout<<"r petsc:\n";
      //  VecView(r,PETSC_VIEWER_STDOUT_WORLD);
      //  std::cout<<"r_input system:\n";
      //  R_input.print(std::cout);
        */
      //  std::cout<<"compute residual done\n"<<std::endl;
        return ierr;
    }

    PetscErrorCode sh_petsc_jac (TS ts, PetscReal t, Vec u, Vec u_t, PetscReal a, Mat *j, Mat *pc, MatStructure *flag, void *ctx)
    {
        int ierr=0;
        //std::cout<<"calling jacobian\n"<<std::endl;

        double tt=static_cast<double>(t);
        double aa=static_cast<double>(a);
        
        LinearImplicitSystem &sys=
        *(static_cast<LinearImplicitSystem*> (ctx));
        
        SH_jacobian(ts,tt,u,u_t,aa,j,pc,ctx);
        
        /*
        std::cout<<"u_t jacob petsc:\n";
        VecView(u_t,PETSC_VIEWER_STDOUT_WORLD);
        
        
        PetscVector<Number> U_input(u), U_t_input(u_t);
        PetscVector<Number>& U_system = *dynamic_cast<PetscVector<Number>*>(sys.solution.get());
        PetscVector<Number>& U_t_system = *dynamic_cast<PetscVector<Number>*>(&(sys.get_vector("u_t")));
        PetscMatrix<Number> J_input(*j);
        PetscMatrix<Number>& J_system = *dynamic_cast<PetscMatrix<Number>*>(sys.matrix);
        
        U_input.swap(U_system);
        U_t_input.swap(U_t_system);
        J_input.swap(J_system);
        
        sys.get_dof_map().enforce_constraints_exactly(sys);
        sys.update();
        SH_jacobian(sys,tt,aa);
        
        // Swap back
        J_system.close();
        U_input.swap(U_system);
        U_t_input.swap(U_t_system);
        J_input.swap(J_system);
        */
       // std::cout<<"compute jacobian done\n"<<std::endl;
        return ierr;
    }
    
    PetscErrorCode sh_petsc_output_monitor(TS ts, PetscInt it_number,PetscReal c_time,Vec u,void *ctx)
    {
        LinearImplicitSystem &sys=
        *(static_cast<LinearImplicitSystem*> (ctx));
        int ploteach=sys.get_equation_systems().parameters.get<int>("ploteach");
        if (it_number%ploteach==0)
        {
            PetscVector<Number> U(u);
            PetscVector<Number>& U_system = *dynamic_cast<PetscVector<Number>*>(sys.solution.get());
            U_system = U;
          //  U_input.swap(U_system);
            sys.update();
            EquationSystems& es = sys.get_equation_systems();
            std::string datadir = es.parameters.get<std::string>("datadir");
            char buffer[50];
            sprintf(buffer, "SHSurf_%09.0f",  c_time);
            OStringStream     filename;
            filename << datadir<<std::string(buffer)<<".pvtu";
            EquationSystems& eq = sys.get_equation_systems();
            VTKIO (eq.get_mesh()).write_equation_systems (filename.str(), eq);
          //  U_input.swap(U_system);
	   
	    PetscViewer petsc_viewer;
	    int ierr = PetscViewerCreate(libMesh::COMM_WORLD, &petsc_viewer); CHKERRABORT(libMesh::COMM_WORLD,ierr);
	//	ierr = PetscViewerSetFormat(petsc_viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRABORT(libMesh::COMM_WORLD,ierr);
	//	ierr=	PetscViewerASCIIOpen(libMesh::COMM_WORLD, "solution.dat", &petsc_viewer); CHKERRABORT(libMesh::COMM_WORLD,ierr);
	    ierr= PetscViewerBinaryOpen(libMesh::COMM_WORLD, "solution.dat", FILE_MODE_WRITE, &petsc_viewer); 
	    CHKERRABORT(libMesh::COMM_WORLD,ierr);
	    ierr = VecView(u,petsc_viewer);
	    PetscViewerDestroy(&petsc_viewer);
	    std::ofstream datafile;
	    datafile.open("solution_data.dat", std::ios::out);
	    datafile<<"time = "<<c_time<<"\n";
	    datafile.close();
        }

    }

    PetscErrorCode sh_StatsMonitor(TS ts, PetscInt it_number,PetscReal c_time,Vec u,void *ctx)
    {
        LinearImplicitSystem &sys=
        *(static_cast<LinearImplicitSystem*> (ctx));
        double free_energy=sys.get_equation_systems().parameters.get<Real>("free_energy");
        PetscReal dt;
        TSGetTimeStep(ts,&dt);
        std::cout<<"Time: "<<c_time<<", dt: "<<dt<<", Free energy U="<<std::setprecision(9)<<free_energy<<std::endl;
        
    }
    
}



// The main program.
int main (int argc, char** argv)
{
    LibMeshInit init (argc, argv);
    GetPot infile("sh.in");
    GetPot command_line(argc, argv);
    bool read_solution=false;
    double starttime=0.0;

    if (command_line.search(1,"-read_solution"))
      {
	read_solution=true;
      }

    std::string input_mesh = infile("Input_Mesh", "input.off");
    double dt = infile("dt", 0.1);
    int Nsteps = infile("nsteps", 1000);

    std::srand(infile("seed", 0));
    
    const unsigned int dim = 2;
    Mesh mesh(dim);
    mesh.read(input_mesh);

	MeshRefinement mesh_refinement (mesh);
    mesh_refinement.uniformly_refine(infile("nrefine",0));
    MeshTools::Modification::flatten(mesh);
	MeshTools::Modification::all_tri(mesh);
    project_to_1_sphere(mesh);
    const Real L = infile("meshscalefactor", 1.0);
    std::cout<<"Scaling mesh by factor L="<<L<<std::endl;
    MeshTools::Modification::scale(mesh, L, L, L);
	MeshTools::Subdiv::prepare_subdiv_mesh(mesh,infile("ghosted",false));
	mesh.print_info();

    EquationSystems equation_systems (mesh);
    
    LinearImplicitSystem & system = equation_systems.add_system<LinearImplicitSystem> ("SH");
    AutoPtr<DofMap> subdivdofmap( new DofMapSubdiv(system.number()));
    system.replace_dof_map(subdivdofmap);
    system.add_variable ("u", FOURTH, SUBDIV);
    system.add_vector ("u_t");
    system.add_vector ("u0");
    
    equation_systems.parameters.set<Real>("a") = infile("a", 0.0);
    equation_systems.parameters.set<Real>("b") = infile("b", 0.0);
    equation_systems.parameters.set<Real>("c") = infile("c", 0.0);
    
    equation_systems.parameters.set<Real>("gamma0") = infile("gamma0", 0.0);
    equation_systems.parameters.set<Real>("free_energy") = 0.0;

    equation_systems.parameters.set<Real>("gamma2") = infile("gamma2", 1.0);
    equation_systems.parameters.set<Real>("penalty") = infile("penalty", 1000);
    equation_systems.parameters.set<std::string>("datadir") = infile("datadir","./");
    equation_systems.parameters.set<int>("ploteach") = infile("ploteach",1);
    
    system.attach_init_function(random_initial_conditions);
    //system.attach_init_function(planar_stripes_initial_conditions);
  //  system.attach_init_function(winded_stripes_initial_conditions);

    equation_systems.init ();
    equation_systems.print_info();
    init_cached_elements(ELEMENT_CACHE, system);

    
    int ierr=0;
    TS ts;
    ierr = TSCreate(libMesh::COMM_WORLD,&ts);CHKERRABORT(libMesh::COMM_WORLD,ierr);
    SparseMatrix<double>&  jac=*(system.matrix);
    NumericVector<double>& sol=*(system.solution);
    NumericVector<double>& res=*(system.rhs);
    
    jac.close(); sol.close(); res.close();
    
    PetscMatrix<double>* J = dynamic_cast<PetscMatrix<double>*>(&jac);
    PetscVector<double>* U   = dynamic_cast<PetscVector<double>*>(&sol);
    PetscVector<double>* R   = dynamic_cast<PetscVector<double>*>(&res);
    
   if (read_solution)
      {
	PetscViewer petsc_viewer;
	ierr = PetscViewerCreate(libMesh::COMM_WORLD, &petsc_viewer); CHKERRABORT(libMesh::COMM_WORLD,ierr);
	//	ierr = PetscViewerSetFormat(petsc_viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRABORT(libMesh::COMM_WORLD,ierr);
	//	ierr=	PetscViewerASCIIOpen(libMesh::COMM_WORLD, "solution.dat", &petsc_viewer); CHKERRABORT(libMesh::COMM_WORLD,ierr);
	ierr= PetscViewerBinaryOpen(libMesh::COMM_WORLD, "solution.dat", FILE_MODE_READ, &petsc_viewer); 
	CHKERRABORT(libMesh::COMM_WORLD,ierr);
	ierr = VecLoad(U->vec(),petsc_viewer);
	PetscViewerDestroy(&petsc_viewer);
	VecView(U->vec(), PETSC_VIEWER_STDOUT_WORLD);
	GetPot datafile("solution_data.dat");
	starttime = datafile("time", 0.0);
      }


    ierr = TSSetSolution(ts,U->vec());CHKERRQ(ierr);
    ierr = TSSetIFunction(ts,R->vec(),sh_petsc_residual,&system);CHKERRQ(ierr);
    ierr = TSSetIJacobian(ts,J->mat(),J->mat(),sh_petsc_jac,&system);CHKERRQ(ierr);
    ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
    ierr = TSSetDuration(ts,Nsteps,1000000.0);CHKERRQ(ierr);
    ierr = TSSetTime(ts,starttime);CHKERRQ(ierr);
    ierr = TSSetTimeStep(ts,dt);CHKERRQ(ierr);
    ierr = TSSetProblemType(ts,TS_LINEAR);CHKERRQ(ierr);
    ierr = TSMonitorSet(ts,sh_petsc_output_monitor,&system,NULL);CHKERRQ(ierr);
    ierr = TSMonitorSet(ts,sh_StatsMonitor,&system,NULL);CHKERRQ(ierr);
    ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
    PetscReal ftime;
    std::cout<<"Solving...\n";

    ierr = TSSolve(ts,U->vec(),&ftime);CHKERRQ(ierr);
    
    ierr = TSDestroy(&ts);CHKERRQ(ierr);

    
    std::cout<<"System solved\n";
	OStringStream file_name,dump_name;
	file_name<<"final-solution.pvtu";
    dump_name<<"solution.dat";
    equation_systems.write (dump_name.str(), libMeshEnums::WRITE);
	VTKIO (mesh).write_equation_systems(file_name.str(), equation_systems);
	
	return 0;

}
