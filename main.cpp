//
// vdwtest.cpp
//
// simple short-range MD test program
//
// Jun Makino March 9 2015
//
// Known problem as of Mar 13 2015
//   -- acc and phi are consistent only for the case of m=1
// This has been fixed as of Mar 15 2015


#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>
#include<particle_mesh.hpp>
#include<param.h>
#include<param_fdps.h>

#include <string>
#include <vector>

#define NATOMTYPE 2
#define NMOLTYPE 2

const PS::F64 epsilon[NATOMTYPE][NATOMTYPE] = {{1.0,1.0},{1.0,1.0}};
const PS::F64 sigma[NATOMTYPE][NATOMTYPE]   = {{1.0,1.0},{1.0,1.0}};

class FileHeader{
public:
  PS::F64vec s,e;

  FileHeader(PS::F64vec _s,PS::F64vec _e) : s(_s),e(_e) {}
  PS::S32 readAscii(FILE * fp){
    return 0;
  }
  void writeAscii(FILE* fp) const{
    fprintf(fp, "'box_sx=%lf,box_sy=%lf,box_sz=%lf,box_ex=%lf,box_ey=%lf,box_ez=%lf\n",
	    s.x,s.y,s.z, e.x,e.y,e.z);
  }
};

class Force{
public:
  PS::F64vec3 acc;
  PS::F64 pot;
  void clear(){
    acc = 0.0;
    pot = 0.0;
  }
};

class FP{
public:
  PS::S64 id;
  PS::S64 type;
  PS::F64 mass;
  PS::F64 charge;
  PS::F64vec3 gpos;
  PS::F64vec3 pos;

  PS::F64vec3 vel;
  PS::F64vec3 acc;

  PS::F64 pot;
  PS::F64 search_radius;

  PS::F64 getRsearch() const {
    return this->search_radius;
  }
  PS::F64vec getPos() const { return gpos; }
  PS::F64 getChargeParticleMesh() const {return charge;};

  void copyFromForce(const Force & force){
    acc = force.acc;
    pot = force.pot;
  }

  void writeAscii(FILE* fp) const{
    fprintf(fp, "%lld %lld %lf %lf %lf\n",
	    id, type, pos.x, pos.y, pos.z);
  }
  void readAscii(FILE* fp){
    fscanf(fp, "%lld %lld %lf %lf %lf\n",
	   &id, &type, &pos.x, &pos.y, &pos.z);
  }

  void IntegratePos(const PS::F64 dt,const PS::F64 cell_size){
    pos += vel*dt/cell_size;
  }
  void IntegrateVel(const PS::F64 dth){
    vel += acc/mass*dth;
  }
};

class EP{
public:
  PS::S64 id;
  PS::S64 type;
  PS::F64vec pos;
  PS::F64 charge;

  PS::F64 search_radius;
  PS::F64 getRSearch() const {
    return this->search_radius;
  }
  PS::F64vec getPos() const { return pos;}
  void copyFromFP(const FP & fp){
    pos = fp.pos;
    id = fp.id;
    type = fp.type;
    charge = fp.charge;
    search_radius = fp.search_radius;
  }
  PS::F64 getCharge() const { return charge; }

  void setPos(const PS::F64vec3 _pos){
    pos = _pos;
  }
};


struct CalcForceEpEp{
  const PS::F64 rcut;
  const PS::F64 cell_size;
  CalcForceEpEp(const PS::F64 _rc,const PS::F64 _cs)
    : rcut(_rc),cell_size(_cs) {}

  void operator () (const EP * ep_i,
		    const PS::S32 n_ip,
		    const EP * ep_j,
		    const PS::S32 n_jp,
		    Force * force){
    const PS::F64 rcsq = rcut*rcut;
    for(PS::S32 i=0; i<n_ip; i++){
      const PS::F64vec ri = ep_i[i].pos;
      PS::F64vec ai = 0.0;
      PS::F64 poti = 0.0;
      const PS::S64 idi = ep_i[i].id;
      const PS::F64 qi = ep_i[i].charge;
      for(PS::S32 j=0; j<n_jp; j++){
	if( idi == ep_j[j].id ) continue;
	const PS::F64vec rij = (ri - ep_j[j].pos)*cell_size;
	const PS::F64 r2 = rij*rij;
	if (r2 < rcsq){
	  const PS::F64 r2_inv = 1.0/r2;
	  const PS::F64 r6_inv = r2_inv * r2_inv * r2_inv;
	  const PS::F64 r12_inv = r6_inv * r6_inv;
	  poti += 4.0*(r12_inv - r6_inv);
	  ai += (48.0*r12_inv - 24.0*r6_inv) * r2_inv * rij;
	  /*
	  const PS::F64 r = sqrt(r2);
	  const PS::F64 r_inv = r2_inv * r;
	  const PS::F64 qqri = qi * ep_j[j].charge * r_inv;
	  poti += qqri;
	  ai += (qqri * r2_inv) * rij;
	  //*/
	}
      }
      force[i].acc += ai;
      force[i].pot += poti;
    }
  }
};

struct CalcForceEpSp{
  CalcForceEpSp(){}

  void operator () (const EP * ep_i,
		    const PS::S32 n_ip,
		    const PS::SPJMonopoleCutoff * ep_j,
		    const PS::S32 n_jp,
		    Force * force){
    for(PS::S32 i=0; i<n_ip; i++){
      const PS::F64vec ri = ep_i[i].pos;
      PS::F64vec ai = 0.0;
      PS::F64 poti = 0.0;
      const PS::F64 qi = ep_i[i].charge;
      for(PS::S32 j=0; j<n_jp; j++){
	const PS::F64vec rij = ri - ep_j[j].pos;
	const PS::F64 r2 = rij*rij;
	const PS::F64 r_inv = 1. / sqrt(r2);

	const PS::F64 qqri = qi * ep_j[j].mass * r_inv;
	poti += qqri;
	ai += (qqri * r_inv * r_inv) * rij;
      }
      force[i].acc += ai;
      force[i].pot += poti;
    }
  }
};


//#define ENABLE_LONG_RANGE
class ForceCalculator {
public:
  //PS::TreeForForceLong<Force, EP, EP>::MonopoleWithCutoff pp;
  PS::TreeForForceShort<Force, EP, EP>::Scatter pp;
#ifdef ENABLE_LONG_RANGE
  PS::PM::ParticleMesh pm;
#endif
  // Methods for water simulation
  template <class Tpsys>
  void initialize(const Tpsys &system) {
    PS::S32 numPtclLocal = system.getNumberOfParticleLocal();
    PS::U64 ntot = 3 * numPtclLocal;
    pp.initialize(ntot,0.0);
  };

  template <class Tpsys,class Tdinfo>
  void operator ()(Tpsys &system,Tdinfo &dinfo,const PS::F64 rcut,const PS::F64 cell_size) {
    //* Local variables
    PS::S32 numPtclLocal = system.getNumberOfParticleLocal();

    //* Reset potential and accelerations
    for (PS::S32 i=0; i<numPtclLocal; i++) {
      system[i].pot  = 0.0;
      system[i].acc = 0.0;
    }
#ifdef ENABLE_LONG_RANGE
    //=================
    //* PM part 
    //=================
    pm.setDomainInfoParticleMesh(dinfo);
    pm.setParticleParticleMesh(system);
    pm.calcMeshForceOnly();
    for (PS::S32 i=0; i<numPtclLocal; i++) { 
      const PS::F32vec pos = system[i].pos; 
      system[i].pot -= pm.getPotential(pos);
      system[i].acc -= pm.getForce(pos);
    }
#endif
    //=================
    //* PP part 
    //=================
    //*
    pp.calcForceAll(CalcForceEpEp(rcut,cell_size),
		    //		    CalcForceEpSp(),
		    system, dinfo);
    //*/
    for (PS::S32 i=0; i<numPtclLocal; i++) {
      const Force result = pp.getForce(i);
      system[i].pot += result.pot;
      system[i].acc += result.acc;
    }
  }
};

void RATTLE_c(FP *atoms,
	      const PS::F64vec3 r_old[3],
	      const PS::F64 cell_size,
	      const PS::F64 dt){
  // iteration parameter
  const int max_iteration = 1000;
  const double max_error = 1e-8;
  // bond parameter
  const int nbond = 3;
  const PS::S64 pair[2*nbond] = {0,1, 1,2, 2,0};
  const PS::F64 bond[nbond] = {0.3,0.3,0.3}; // test

  const double dti = 1.0/dt;
  for(int iter=1;iter<=max_iteration;iter++){
    bool isEnough = true;
    for(int b=0;b<nbond;b++){
      const int i = pair[2*b+0];
      const int j = pair[2*b+1];
      const PS::F64vec3 rij = (atoms[i].pos - atoms[j].pos)*cell_size;
      const PS::F64 r2 = rij*rij;
      const PS::F64 b2 = bond[b]*bond[b];
      if(fabs((sqrt(r2) - bond[b])/bond[b]) > max_error){
	isEnough = false;
	const PS::F64 mi_i = 1.0 / atoms[i].mass;
	const PS::F64 mi_j = 1.0 / atoms[j].mass;
	const PS::F64vec3 rij_old = (r_old[i] - r_old[j])*cell_size;
	const PS::F64 gij = (r2 - b2) / (2.0*(rij*rij_old) * (mi_i + mi_j));
	atoms[i].pos -= (gij*mi_i * rij_old)/cell_size;
	atoms[j].pos += (gij*mi_j * rij_old)/cell_size;
	atoms[i].vel -= gij*mi_i*dti * rij_old;
	atoms[j].vel += gij*mi_j*dti * rij_old;
      }
    }
    if(isEnough) break;
    if(iter==max_iteration){
      std::cerr << "warning: SHAKE iteration is too long (>" << max_iteration << ")" << std::endl;
    }
  }
}

void RATTLE_v(FP *atoms,
	      const PS::F64 cell_size,
	      const PS::F64 dt){
  // iteration parameter
  const int max_iteration = 1000;
  const double max_error = 1e-12;
  // bond parameter
  const int nbond = 3;
  const PS::S64 pair[2*nbond] = {0,1, 1,2, 2,0};
  const PS::F64 bond[nbond] = {0.3,0.3,0.3}; // test

  for(int iter=1;iter<=max_iteration;iter++){
    bool isEnough = true;
    for(int b=0;b<nbond;b++){
      const int i = pair[2*b+0];
      const int j = pair[2*b+1];
      const PS::F64vec3 rij = (atoms[i].pos - atoms[j].pos)*cell_size;
      const PS::F64vec3 vij = atoms[i].vel - atoms[j].vel;
      const PS::F64 rv = rij*vij;
      const PS::F64 delta = rv*dt/bond[b];
      if(fabs(delta) > max_error){
	isEnough = false;
	const PS::F64 b2 = bond[b]*bond[b];
	const PS::F64 mi_i = 1.0 / atoms[i].mass;
	const PS::F64 mi_j = 1.0 / atoms[j].mass;
	const PS::F64 kij = rv / (b2 * (mi_i + mi_j));
	atoms[i].vel -= kij*mi_i * rij;
	atoms[j].vel += kij*mi_j * rij;
      }
    }
    if(isEnough) break;
    if(iter==max_iteration){
      std::cerr << "warning: RATTLE iteration is too long (>" << max_iteration << ")" << std::endl;
    }
  }
}

void MakeFaceCubicCenter(const long long int n_tot,
			 double *&mass,
			 PS::F64vec *&pos,
			 PS::F64vec *&vel,
			 const double density,
			 const double eng = -0.25,
			 const int seed = 0){
  assert(eng < 0.0);
  //static const double PI = atan(1.0) * 4.0;
  mass = new double[n_tot];
  pos = new PS::F64vec3[n_tot];
  vel = new PS::F64vec3[n_tot];

  PS::MTTS mt;
  const int nmol = n_tot / 3;
  int nunit = 1;
  while(4*nunit*nunit*nunit < nmol) nunit++;
  if (nmol != 4*nunit*nunit*nunit){
    std::cerr << "MakeFaceCubicCenter: n_tot and 4*nunit^3 must be the same. "
	      << n_tot << "!= " << 4*nunit*nunit*nunit <<std::endl;
    PS::Abort();
  }
  mt.init_genrand(PS::Comm::getRank()*PS::Comm::getNumberOfThread()+PS::Comm::getThreadNum());

  PS::F64 unit_size = 1.0/(PS::F64)nunit;
  printf("unit_size is %lf\n",unit_size);
  PS::F64 ush = unit_size * 0.5;
  PS::F64vec unit[4];
  unit[0].x = 0.0; unit[1].x = ush; unit[2].x = 0.0; unit[3].x = ush;
  unit[0].y = 0.0; unit[1].y = ush; unit[2].y = ush; unit[3].y = 0.0;
  unit[0].z = 0.0; unit[1].z = 0.0; unit[2].z = ush; unit[3].z = ush;

  int ip=0;
  for(int i=0; i<nunit; i++){
    for(int j=0; j<nunit; j++){
      for(int k=0; k<nunit; k++){
	for(int l=0; l<4; l++){
	  pos[3*ip+0].x = i*unit_size + unit[l].x;
	  pos[3*ip+0].y = j*unit_size + unit[l].y;
	  pos[3*ip+0].z = k*unit_size + unit[l].z;

	  pos[3*ip+1].x = pos[3*ip+0].x + ush*0.1;
	  pos[3*ip+1].y = pos[3*ip+0].y;
	  pos[3*ip+1].z = pos[3*ip+0].z;

	  pos[3*ip+2].x = pos[3*ip+0].x;
	  pos[3*ip+2].y = pos[3*ip+0].y + ush*0.1;
	  pos[3*ip+2].z = pos[3*ip+0].z;
	  ip++;
	}
      }
    }
  }
  assert(ip == nmol && 3*ip == n_tot);

  for(int i=0; i<nmol; i++){
    mass[3*i+0] = mass[3*i+1] = mass[3*i+2] = 1.0;
    const double v_max = 0.1;
    do {
      vel[3*i+0][0] = (2. * mt.genrand_res53() - 1.) * v_max;
      vel[3*i+0][1] = (2. * mt.genrand_res53() - 1.) * v_max;
      vel[3*i+0][2] = (2. * mt.genrand_res53() - 1.) * v_max;

      vel[3*i+1][0] = vel[3*i][0];
      vel[3*i+1][1] = vel[3*i][1];
      vel[3*i+1][2] = vel[3*i][2];

      vel[3*i+2][0] = vel[3*i][0];
      vel[3*i+2][1] = vel[3*i][1];
      vel[3*i+2][2] = vel[3*i][2];
    }while(vel[3*i] * vel[3*i] >= v_max * v_max);
  }

  PS::F64vec cm_vel = 0.0;
  double  cm_mass = 0.0;
  for(int i=0; i<n_tot; i++){
    cm_vel += mass[i] * vel[i];
    cm_mass += mass[i];
  }
  cm_vel /= cm_mass;
  for(int i=0; i<n_tot; i++){
    vel[i] -= cm_vel;
  }
}

template<class Tpsys>
void SetParticlesFCC(Tpsys & psys,
		     const PS::S32 n_tot,
		     const double density,
		     const double temperature){
  PS::F64 * mass;
  PS::F64vec * pos;
  PS::F64vec * vel;

  const PS::F64 eng = -0.25;
  const int nmol = n_tot/3;
  MakeFaceCubicCenter(n_tot, mass, pos, vel, density, eng);
  PS::F64vec *gpos = new PS::F64vec[nmol];
  for(int i=0; i<nmol; i++){
    gpos[i] = 0.0;
    PS::F64 m = 0.0;
    for(int j=0;j<3;j++){
      gpos[i] += mass[3*i+j]*pos[3*i+j];
      m += mass[3*i+j];
    }
    gpos[i] /= m;
    //printf("gpos %d: %lf %lf %lf\n",i,gpos[i].x,gpos[i].y,gpos[i].z);
  }

  PS::S32 n_proc = PS::Comm::getNumberOfProc();
  PS::S32 rank = PS::Comm::getRank();
  printf("rank %d: n_proc = %d\n",rank,n_proc);

  PS::S32 nmol_loc = nmol / n_proc;
  if(nmol%n_proc > rank) nmol_loc++;
  PS::S32 i_h = (nmol/n_proc)*rank;
  if(nmol%n_proc > rank) i_h += rank;

  psys.setNumberOfParticleLocal(3*nmol_loc);
  const PS::F64 cell_size = powf((PS::F64)nmol/density,1.0/3.0);
  for(int i=0; i<nmol_loc; i++){
    for(int j=0;j<3;j++){
      const int mid = i_h+i;
      const int aid = 3*mid + j;
      psys[3*i+j].mass = mass[aid];
      psys[3*i+j].gpos = gpos[mid];
      psys[3*i+j].pos = pos[aid];
      psys[3*i+j].vel = vel[aid];
      psys[3*i+j].id = mid;
      psys[3*i+j].type = j%3==0?0:1;
      psys[3*i+j].search_radius = 5.0/cell_size;
      assert(0.0<=psys[3*i+j].gpos.x && psys[3*i+j].gpos.x<1.0);
      assert(0.0<=psys[3*i+j].gpos.y && psys[3*i+j].gpos.y<1.0);
      assert(0.0<=psys[3*i+j].gpos.z && psys[3*i+j].gpos.z<1.0);
    }
  }
  delete [] mass;
  delete [] pos;
  delete [] vel;

  ScaleVelocity(psys,temperature);
}

//#define TEST
template<class Tpsys>
void SetGravityCenterForAtoms(Tpsys & system){
  const PS::S32 n = system.getNumberOfParticleLocal();
  const unsigned int nmol = n/3;
  for(unsigned int i=0;i<nmol;i++){
    PS::F64vec3 gpos = 0.0;
    PS::F64 mass = 0.0;
    for(int j=0;j<3;j++){
      gpos += system[3*i+j].mass * system[3*i+j].pos;
      mass += system[3*i+j].mass;
    }
    gpos /= mass;
    for(int j=0;j<3;j++) system[3*i+j].gpos = gpos;
  }
}

template<class Tpsys>
void PeriodicBoundaryCondition(Tpsys & system,
			       const PS::F64 boxdh){
  PS::S32 n = system.getNumberOfParticleLocal();
  for(int i=0; i<n; i++){
    for(int k=0;k<3;k++){
      if(system[i].gpos[k] <  0.0){
	system[i].pos[k]  += 1.0;
	system[i].gpos[k] += 1.0;
      }
      if(system[i].gpos[k] >= 1.0){
	system[i].pos[k]  -= 1.0;
	system[i].gpos[k] -= 1.0;
      }
    }
  }
}

template<class Tpsys>
void RemoveTotalMomentum(Tpsys &system){
  const PS::S32 n_loc = system.getNumberOfParticleLocal();
  PS::F64vec cm_vel_loc = 0.0;
  PS::F64    cm_mass_loc = 0.0;
  for(int i=0; i<n_loc; i++){
    cm_vel_loc += system[i].mass * system[i].vel;
    cm_mass_loc += system[i].mass;
  }
  PS::F64vec cm_vel;
  PS::F64    cm_mass;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&cm_vel_loc.x, &cm_vel.x, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&cm_vel_loc.y, &cm_vel.y, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&cm_vel_loc.z, &cm_vel.z, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&cm_mass_loc,  &cm_mass,  1, PS::GetDataType<PS::F64>(), MPI::SUM);
#else
  cm_vel = cm_vel_loc;
  cm_mass = cm_mass_loc;
#endif
  cm_vel /= cm_mass;
  for(int i=0; i<n_loc; i++){
    system[i].vel -= cm_vel;
  }
}

template<class Tpsys>
void ScaleVelocity(Tpsys & system,const PS::F64 T){
  const PS::S32 natom_local = system.getNumberOfParticleLocal();
  PS::F64 ekin_loc = 0.0;
  for(PS::S32 i=0; i<natom_local; i++){
    ekin_loc += system[i].mass * system[i].vel * system[i].vel;
  }
  ekin_loc *= 0.5;
  PS::S32 natom;
  PS::F64 ekin;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&ekin_loc, &ekin, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&natom_local, &natom, 1, PS::GetDataType<PS::S32>(), MPI::SUM);
#else
  ekin = ekin_loc;
  natom = natom_local;
#endif
  const PS::F64 scaler = sqrt(1.5*natom*T / ekin);
  for(PS::S32 i=0;i<natom_local;i++) system[i].vel *= scaler;

  RemoveTotalMomentum(system);
}

template<class Tpsys>
void CalcEnergy(const Tpsys & system,
                PS::F64 & etot,
                PS::F64 & ekin,
                PS::F64 & epot,
                const bool clear=true){
  if(clear){
    etot = ekin = epot = 0.0;
  }
  PS::F64 etot_loc = 0.0;
  PS::F64 ekin_loc = 0.0;
  PS::F64 epot_loc = 0.0;
  const PS::S32 nbody = system.getNumberOfParticleLocal();
  for(PS::S32 i=0; i<nbody; i++){
    ekin_loc += system[i].mass * system[i].vel * system[i].vel;
    epot_loc += system[i].pot;
  }
  ekin_loc *= 0.5;
  epot_loc *= 0.5;
  etot_loc = ekin_loc + epot_loc;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&etot_loc, &etot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&epot_loc, &epot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&ekin_loc, &ekin, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
#else
  etot = etot_loc;
  epot = epot_loc;
  ekin = ekin_loc;
#endif
}

int main(int argc, char *argv[]){
  PS::Initialize(argc, argv);

  PS::F64 Tbegin = PS::GetWtime();
  std::cout<<std::setprecision(15);
  std::cerr<<std::setprecision(15);
  PS::F64 theta = 0.5;
  const PS::S32 n_leaf_limit = 8;

  long long int nmol = 500;
  PS::F64 density     = 1.05;
  PS::F64 temperature = 0.8;

  PS::F64 dt = 0.0001;
  PS::F64 rcut = 4.5;

  PS::S32 n_group_limit = 64;
  PS::S32 nstep     = 1000;
  PS::S32 nstep_eq  = 1000;
  PS::S32 nstep_snp = 100;
  PS::S32 nstep_diag = 100;

  //char sinput[1024];
  char dir_name[1024];
  int c;
  sprintf(dir_name,"./result");
  while((c=getopt(argc,argv,"o:N:d:T:s:e:S:D:t:c:n:h")) != -1){
    switch(c){
    case 'o':
      sprintf(dir_name,optarg);
      break;
    case 'N':
      nmol = atoi(optarg);
      std::cerr<<"nmol="<<nmol<<std::endl;
      break;
    case 'd':
      density = atof(optarg);
      std::cerr<<"density="<<density<<std::endl;
      break;
    case 'T':
      temperature = atof(optarg);
      std::cerr<<"temperature="<<temperature<<std::endl;
      break;
    case 's':
      nstep = atoi(optarg);
      std::cerr<<"nstep="<<nstep<<std::endl;
      break;
    case 'e':
      nstep_eq = atoi(optarg);
      std::cerr<<"nstep_eq="<<nstep_eq<<std::endl;
      break;
    case 'S':
      nstep_snp = atoi(optarg);
      std::cerr<<"nstep_snp="<<nstep_snp<<std::endl;
      break;
    case 'D':
      nstep_diag = atoi(optarg);
      std::cerr<<"nstep_diag="<<nstep_diag<<std::endl;
      break;
    case 't':
      dt = atof(optarg);
      std::cerr<<"dt="<<dt<<std::endl;
      break;
    case 'c':
      rcut = atof(optarg);
      std::cerr<<"rcut="<<rcut<<std::endl;
      break;
    case 'n':
      n_group_limit = atoi(optarg);
      std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
      break;
    case 'h':
      std::cerr<<"N: n_tot (default: 1000)"<<std::endl;
      std::cerr<<"d: number density (default: 1.05)"<<std::endl;
      std::cerr<<"T: temperature (default: 0.8)"<<std::endl;
      std::cerr<<"s: number of steps (default: 1000)"<<std::endl;
      std::cerr<<"S: time step for snapshot(default: 100)"<<std::endl;
      std::cerr<<"D: time step for diag(default: 100)"<<std::endl;
      std::cerr<<"e: number of steps for equilibration(default: 1000)"<<std::endl;
      std::cerr<<"o: dir name of output (default: ./result)"<<std::endl;
      std::cerr<<"t: time step (default: 0.0001)"<<std::endl;
      std::cerr<<"n: n_group_limit (default: 64.0)"<<std::endl;
      return 0;
    }
  }
  long long int n_tot = 3*nmol;
  PS::F64 boxdh = 0.5*powf(nmol/density,1./3.);
  if(boxdh < rcut){
    rcut = boxdh;
    std::cerr << "warning: rcut > 0.5*box_size. rcut is set to " << rcut << std::endl;
  }
  const PS::F64 box_size = 2.0 * boxdh;

  struct stat st;
  if(stat(dir_name, &st) != 0) {
    PS::S32 rank = PS::Comm::getRank();
    PS::S32 ret_loc, ret=0;
    if(rank == 0)
      ret_loc = mkdir(dir_name, 0777);
    PS::Comm::broadcast(&ret_loc, ret);
    if(ret == 0) {
      if(rank == 0)
	fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
    } else {
      fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
      PS::Abort();
      exit(0);
    }
  }

  std::ofstream fout_eng;
  std::ofstream fout_tcal;
  if(PS::Comm::getRank() == 0){
    char sout_de[1024];
    char sout_tcal[1024];
    sprintf(sout_de, "%s/t-de.dat", dir_name);
    sprintf(sout_tcal, "%s/t-tcal.dat", dir_name);
    std::cerr<<sout_de<<std::endl;
    std::cerr<<sout_tcal<<std::endl;
    fout_eng.open(sout_de);
    fout_tcal.open(sout_tcal);
  }

  if(PS::Comm::getRank()==0) printf("initializing particle system ...\n");
  PS::ParticleSystem<FP> system_water;
  system_water.initialize();
  if(PS::Comm::getRank()==0) printf("particle system is initialized!\n");

  PS::S32 n_grav_glb = n_tot;
  SetParticlesFCC(system_water, n_tot, density, temperature);
  for(int i=0;i<system_water.getNumberOfParticleLocal()/3;i++){
    PS::F64vec r_old[3] = {system_water[3*i+0].pos,system_water[3*i+1].pos,system_water[3*i+2].pos};
    RATTLE_c(&system_water[3*i],r_old,box_size,dt);
  }
  PeriodicBoundaryCondition(system_water,boxdh);
  if(PS::Comm::getRank()==0) printf("FCC is generated!\n");
  {
    FileHeader header(PS::F64vec(0.0,0.0,0.0),PS::F64vec(1.0,1.0,1.0));
    char filename[256];
    sprintf(filename, "%s/test.cdv", dir_name);
    system_water.writeParticleAscii(filename, header);
  }

  const PS::F64 coef_ema = 0.3;
  PS::DomainInfo dinfo;
  dinfo.initialize(coef_ema);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0,0.0,0.0),
			 PS::F64vec(1.0,1.0,1.0));
  dinfo.collectSampleParticle(system_water);
  dinfo.decomposeDomain();
  if(PS::Comm::getRank()==0) printf("Domain was decomposed\n");
  system_water.exchangeParticle(dinfo);
  if(PS::Comm::getRank()==0) printf("Particles were exchanged\n");

  ForceCalculator pppm;
  pppm.initialize(system_water);
  pppm(system_water,dinfo,rcut,box_size);
  if(PS::Comm::getRank()==0) printf("Force calculation is done\n");

  PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
  ScaleVelocity(system_water,temperature);
  CalcEnergy(system_water, Etot0, Ekin0, Epot0);
  if(PS::Comm::getRank() == 0) printf("Etot = %lf, Epot = %lf, Ekin = %lf\n",Etot0,Epot0,Ekin0);
  PS::F64 Tloop = 0.0;

  PS::S32 snp_id = 0;
  PS::S32 time_snp = 0;
  PS::S32 time_diag = 0;
  bool isInitialized = false;
  PS::F64 time_sys = -dt * nstep_eq;

  PS::F64 Epot_ave = 0.0, Ekin_ave = 0.0;
  int n_loc = system_water.getNumberOfParticleLocal();
  assert(n_loc%3 == 0);
  int nmol_loc = n_loc / 3;
  if(PS::Comm::getRank()==0) printf("Main loop start!\n");
  for(int s=-nstep_eq;s<nstep;s++){

    if(s < 0) ScaleVelocity(system_water,temperature);
    if(s%1000==0) RemoveTotalMomentum(system_water);
    PS::Timer timer;
    timer.reset();
    timer.start();
    if(s == time_snp){
      FileHeader header(PS::F64vec(-boxdh,-boxdh,-boxdh),PS::F64vec(boxdh,boxdh,boxdh));
      char filename[256];
      sprintf(filename, "%s/%04d.cdv", dir_name, snp_id++);
      system_water.writeParticleAscii(filename, header);

      time_snp += nstep_snp;
    }
    timer.restart("WriteParticleAscii");
    if(!isInitialized && s == 0){
      CalcEnergy(system_water, Etot0, Ekin0, Epot0);
      if(PS::Comm::getRank() == 0) printf("Etot0 = %lf, Epot0 = %lf, Ekin0 = %lf\n",Etot0,Epot0,Ekin0);
      isInitialized = true;
    }
    for(int i=0;i<nmol_loc;i++){
      for(int j=0;j<3;j++) system_water[3*i+j].IntegrateVel(0.5*dt);
      PS::F64vec r_old[3] = {system_water[3*i+0].pos,system_water[3*i+1].pos,system_water[3*i+2].pos};
      for(int j=0;j<3;j++) system_water[3*i+j].IntegratePos(dt,box_size);
      RATTLE_c(&system_water[3*i],r_old,box_size,dt);
    }
    timer.restart("SHAKE");
    SetGravityCenterForAtoms(system_water);
    PeriodicBoundaryCondition(system_water,boxdh);
    if(s%10 == 0){
      dinfo.collectSampleParticle(system_water);
      timer.restart("collect");
      dinfo.decomposeDomain();
    }
    else{
      timer.restart("collect");
    }

    timer.restart("decompose");
    system_water.exchangeParticle(dinfo);
    n_loc = system_water.getNumberOfParticleLocal();
    assert(n_loc%3 == 0);
    nmol_loc = n_loc/3;
    Tloop = PS::GetWtime();
    pppm(system_water, dinfo,rcut,box_size);
    Tloop = PS::GetWtime() - Tloop;
	
    for(int i=0;i<nmol_loc;i++){
      for(int j=0;j<3;j++)
	system_water[3*i+j].IntegrateVel(0.5*dt);
      RATTLE_v(&system_water[3*i],box_size,dt);
    }
    timer.stop("Kick");

    fout_tcal<<"time_sys= "<<time_sys<<std::endl;
    fout_tcal<<"tree_water.getMemSizeUsed()= "<<pppm.pp.getMemSizeUsed()*1e-9<<" [Gbyte]";
    fout_tcal<<" system_water.getMemSizeUsed()= "<<system_water.getMemSizeUsed()*1e-9<<" [Gbyte]"<<std::endl;
    fout_tcal<<"Tloop= "<<Tloop<<" Ttot="<<PS::GetWtime()-Tbegin<<std::endl;
    timer.dump(fout_tcal);
    fout_tcal<<std::endl;

    CalcEnergy(system_water, Etot1, Ekin1, Epot1);
    if(s>=0){
      Epot_ave += Epot1;
      Ekin_ave += Ekin1;
    }
    if(s == time_diag) {
      if(PS::Comm::getRank() == 0){
	/*
	ref_grav(system_water,n_tot,boxdh,Epot_ref,Ekin_ref);
	printf("diff of FDPS and ref is %e %e \n",(Epot1-Epot_ref)/Epot1,(Ekin1-Ekin_ref)/Ekin1);
	//*/
	fout_eng<<time_sys<<"   "<< " " << Epot1 << " " << Ekin1 << " " <<(Etot1-Etot0)/Etot0<<std::endl;
	fprintf(stderr, "%10.7f %lf %lf %+e\n",
		time_sys, Epot1, Ekin1, (Etot1 - Etot0) / Etot0);
	time_diag += nstep_diag;
      }
      /*
      for(int i=0;i<15;i++){
	printf("%d: %lld %lf %lf %lf\n",i,system_water[i].id,system_water[i].pos.x,system_water[i].pos.y,system_water[i].pos.z);
      }
      //*/
    }
    time_sys += dt;
  }
  PS::Finalize();
  return 0;
}
