#ifndef _TESTMULTIGRID_H_
#define _TESTMULTIGRID_H_

#include <jsoncpp/json/json.h>
#include <mpi/mpi.h>
#include <fstream>
#include <string>
#include "Multigrid/MGFactory.h"

template <int Dim>
class TestMultigrid;

template <int Dim>
class TestMultigrid {
 public:
  TestMultigrid(const std::string& jsonfile, bool useMpi = false);

  TestMultigrid(const std::string amultigridID,
                const RectDomain<Dim> aDomain,
                const char* aBCTypes,
                bool useMpi = false);

  void setParam(const MGParam& aParam) const;

  Tensor<Real, Dim> solve(const ScalarFunction<Dim>* pfunc,
                          const std::array<ScalarFunction<Dim>*, 3> pfuncs,
                          bool useFMVCycle = false,
                          bool useMpi = false) const;

  Real computeError(const Tensor<Real, Dim>& res,
                    const ScalarFunction<Dim>* pfunc,
                    const int p = 0) const;

  void plot(const Tensor<Real, Dim>& res, const std::string& file) const;

  void test(const ScalarFunction<Dim>* pfunc,
            const std::array<ScalarFunction<Dim>*, 3> pfuncs,
            const int numEncryp,
            bool useFMVCycle = false,
            bool useMpi = false) const;

 protected:
  std::unique_ptr<MultigridSolver<Dim> > pM;
  std::string multigridID;
};

template <int Dim>
TestMultigrid<Dim>::TestMultigrid(const std::string amultigridID,
                                  const RectDomain<Dim> aDomain,
                                  const char* aBCTypes,
                                  bool useMpi) {
  MGFactory<Dim>& MGF = MGFactory<Dim>::createFactory();
  registerFactory<Dim>();
  if (!useMpi) {
    pM = MGF.createMultigridSolver(amultigridID, aDomain, aBCTypes);
  } else {
    int mpiSize;
    int mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    int nx = std::round(std::pow(mpiSize, 1.0 / Dim));
    assert(std::pow(nx, Dim) == mpiSize &&
           "mpi parallel should be Dim-th pow.");
    Vec<int, Dim> cores(nx);
    auto stride = Tensor<int, Dim>::strideFromSize(cores);
    auto coord = Tensor<int, Dim>::index2coord(mpiRank, stride);
    Vec<int, Dim> dist = ((aDomain.size() - 1) / nx);
    Vec<int, Dim> lo = coord * dist;
    Vec<int, Dim> hi = min(aDomain.size() - 1, (coord + 1) * dist);
    Box<Dim> bx(lo, hi);
    Vec<Real, Dim> dx = aDomain.spacing();
    RectDomain<Dim> Domain(bx, dx, NodeCentered, 1);
    char BCTypes[4] = {'G', 'G', 'G', 'G'};
    switch (Dim) {
      case 2:
        if (lo[1] == 0)
          BCTypes[0] = aBCTypes[0];
        if (hi[1] == aDomain.size()[1] - 1)
          BCTypes[1] = aBCTypes[1];
      case 1:
        if (lo[0] == 0)
          BCTypes[2] = aBCTypes[2];
        if (hi[0] == aDomain.size()[0] - 1)
          BCTypes[3] = aBCTypes[3];
    }
    pM = MGF.createMultigridSolver(amultigridID, Domain, BCTypes);
  }
}

template <int Dim>
TestMultigrid<Dim>::TestMultigrid(const std::string& jsonfile, bool useMpi) {
  Json::Value root;
  Json::Reader reader;
  std::ifstream os;
  os.open(jsonfile);
  reader.parse(os, root);
  os.close();
  Json::Value Domain = root["Domain"];
  Vec<Real, Dim> rlo;
  Vec<Real, Dim> rhi;
  switch (Dim) {
    case 1:
      rlo = Vec<Real, Dim>{Domain["lo"][0].asDouble()};
      rhi = Vec<Real, Dim>{Domain["hi"][0].asDouble()};
      break;
    case 2:
      rlo = Vec<Real, Dim>{Domain["lo"][0].asDouble(),
                           Domain["lo"][1].asDouble()};
      rhi = Vec<Real, Dim>{Domain["hi"][0].asDouble(),
                           Domain["hi"][1].asDouble()};
      break;
  }
  int an = Domain["n"].asInt();
  const char* bctypes = root["bctypes"].asString().c_str();
  char BCTypes[4];
  switch (Dim) {
    case 1:
      BCTypes[0] = bctypes[0];
      BCTypes[1] = bctypes[1];
      break;
    case 2:
      BCTypes[0] = bctypes[0];
      BCTypes[1] = bctypes[1];
      BCTypes[2] = bctypes[2];
      BCTypes[3] = bctypes[3];
      break;
  }
  std::string restriction = root["restriction"].asString();
  std::string interpolation = root["interpolation"].asString();
  multigridID = restriction + "+" + interpolation;
  const int numpreiter = root["numiter1"].asInt();
  const int numpostiter = root["numiter2"].asInt();
  const int numbottomiter = root["numiter3"].asInt();
  const int maxiter = root["maxiter"].asInt();
  const Real reltol = root["reltol"].asDouble();
  Vec<int, Dim> lo(0);
  Vec<int, Dim> hi(an);
  Box<Dim> bx(lo, hi);
  Vec<Real, Dim> dx = (rhi - rlo) / an;
  RectDomain<Dim> aDomain(bx, dx, NodeCentered, 1);
  TestMultigrid<Dim> tmp(multigridID, aDomain, BCTypes, useMpi);
  pM = std::move(tmp.pM);
  MGParam aParam{numpreiter, numpostiter, numbottomiter, reltol, maxiter};
  pM->setParam(aParam);
}

template <int Dim>
void TestMultigrid<Dim>::setParam(const MGParam& aParam) const {
  pM->setParam(aParam);
}

template <int Dim>
Tensor<Real, Dim> TestMultigrid<Dim>::solve(
    const ScalarFunction<Dim>* pfunc,
    const std::array<ScalarFunction<Dim>*, 3> pfuncs,
    bool useFMVCycle,
    bool useMpi) const {
  const RectDomain<Dim>& Domain = (pM->getvDomain())[0];
  FuncFiller<Dim> FuncF(Domain);
  GhostFiller<Dim> GhostF(Domain);
  Tensor<Real, Dim> rhs(Domain.getGhostedBox());
  rhs = 0;
  FuncF.fill(rhs, pfunc);
  const char* BCTypes = pM->getBCTypes();
  GhostF.transformInitMpi();
  GhostF.fillAllSidesMpi(rhs, BCTypes, pfuncs);
  GhostF.transformFreeMpi();
  Tensor<Real, Dim> phi(Domain.getGhostedBox());
  phi = 0.0;
  pM->solve(phi, rhs, pfuncs, useFMVCycle, useMpi);
  return phi;
}

template <int Dim>
Real TestMultigrid<Dim>::computeError(const Tensor<Real, Dim>& res,
                                      const ScalarFunction<Dim>* pfunc,
                                      const int p) const {
  const RectDomain<Dim>& Domain = (pM->getvDomain())[0];
  Tensor<Real, Dim> exactres(Domain);
  const Vec<int, Dim>& lo = Domain.lo();
  const Vec<Real, Dim>& dx = Domain.spacing();
  switch (Dim) {
    case 1:
      loop_box_1(Domain, i) {
        Vec<int, Dim> Node{i};
        Vec<Real, Dim> rNode = (Node) * dx;
        exactres(i) = (*pfunc)(rNode);
      }
      break;
    case 2:
      loop_box_2(Domain, i, j) {
        Vec<int, Dim> Node{i, j};
        Vec<Real, Dim> rNode = (Node) * dx;
        exactres(i, j) = (*pfunc)(rNode);
      }
      break;
  }
  Tensor<Real, Dim> Rres = res.slice(Domain);
  Laplacian<Dim> POP(Domain);
  return POP.computeNorm(Rres - exactres, p);
}

template <int Dim>
void TestMultigrid<Dim>::plot(const Tensor<Real, Dim>& res,
                              const std::string& file) const {
  const RectDomain<Dim>& Domain = (pM->getvDomain())[0];
  Tensor<Real, Dim> Rres = res.slice(Domain);
  const Vec<int, Dim>& lo = Domain.lo();
  const Vec<Real, Dim>& dx = Domain.spacing();
  std::ofstream os;
  os.open(file);
  os << "res=[\n" << std::endl;
  switch (Dim) {
    case 1:
      loop_box_1(Domain, i) {
        Vec<int, Dim> Node{i};
        Vec<Real, Dim> rNode = (Node - lo) * dx;
        os << rNode[0] << "," << Rres(i) << ";\n" << std::endl;
      }
      os << "];\n" << std::endl;
      os << "x=res(:,1);y=res(:,2);\n" << std::endl;
      os << "plot(x,y);\n" << std::endl;
      break;
    case 2:
      loop_box_2(Domain, i, j) {
        Vec<int, Dim> Node{i, j};
        Vec<Real, Dim> rNode = (Node - lo) * dx;
        os << rNode[0] << "," << rNode[1] << "," << Rres(i, j) << ";\n"
           << std::endl;
      }
      os << "];\n" << std::endl;
      os << "x=res(:,1);y=res(:,2);z=res(:,3);\n" << std::endl;
      os << "T = delaunay(x,y);trimesh(T,x,y,z);\n" << std::endl;
      break;
  }
  os.close();
}

template <int Dim>
void TestMultigrid<Dim>::test(const ScalarFunction<Dim>* pfunc,
                              const std::array<ScalarFunction<Dim>*, 3> pfuncs,
                              const int numEncryp,
                              bool useFMVCycle,
                              bool useMpi) const {
  std::vector<Real> res1(2 * numEncryp + 1);
  std::vector<Real> res2(2 * numEncryp + 1);
  std::vector<Real> res3(2 * numEncryp + 1);
  const char* BCTypes = pM->getBCTypes();
  const RectDomain<Dim>& Domain = (pM->getvDomain())[0];
  int n = Domain.size()[0] - 1;
  const MGParam& Param = pM->getMGParam();
  int rank;
  if (useMpi){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    pM->fillerInitMpi();
  }
  if (!useMpi || rank == 0)
    switch (Dim) {
      case 1:
        std::cout << "Grid " << n << ":" << std::endl;
        break;
      case 2:
        std::cout << "Grid " << n << "X" << n << ":" << std::endl;
        break;
    }
  Tensor<Real, Dim> phi = this->solve(pfunc, pfuncs, useFMVCycle, useMpi);
  res1[0] = this->computeError(phi, pfuncs[0], 0);
  res2[0] = this->computeError(phi, pfuncs[0], 1);
  res3[0] = this->computeError(phi, pfuncs[0], 2);
  std::vector<RectDomain<Dim> > vDomain{Domain};
  for (int i = 1; i <= numEncryp; i++) {
    RectDomain<Dim> nDomain = vDomain.back().refine();
    vDomain.push_back(nDomain);
    TestMultigrid<Dim> nTM(multigridID, nDomain, BCTypes);
    if (useMpi) {
      nTM.pM->fillerInitMpi();
    }
    nTM.setParam(Param);
    n *= 2;
    if (!useMpi || rank == 0)
      switch (Dim) {
        case 1:
          std::cout << "Grid " << n << ":" << std::endl;
          break;
        case 2:
          std::cout << "Grid " << n << "X" << n << ":" << std::endl;
          break;
      }
    Tensor<Real, Dim> cphi = nTM.solve(pfunc, pfuncs, useFMVCycle, useMpi);
    res1[2 * i] = nTM.computeError(cphi, pfuncs[0], 0);
    res2[2 * i] = nTM.computeError(cphi, pfuncs[0], 1);
    res3[2 * i] = nTM.computeError(cphi, pfuncs[0], 2);
    res1[2 * i - 1] = log(res1[2 * i - 2] / res1[2 * i]) / log(2);
    res2[2 * i - 1] = log(res2[2 * i - 2] / res2[2 * i]) / log(2);
    res3[2 * i - 1] = log(res3[2 * i - 2] / res3[2 * i]) / log(2);
    if (useMpi) {
      nTM.pM->fillerFreeMpi();
    }
  }
   if (useMpi){
    pM->fillerFreeMpi();
  }
  if (!useMpi) {
    std::cout << "----------------------------------------" << std::endl;
    const Vec<Real, Dim>& dx = Domain.spacing();
    const Vec<int, Dim>& lo = Domain.lo();
    const Vec<int, Dim>& hi = Domain.hi();
    n = hi[0] - lo[0];
    switch (Dim) {
      case 1:
        std::cout << "\\noindent Domain: "
                  << "$(" << lo[0] << "," << (hi[0] - lo[0]) * dx[0] << ")$"
                  << std::endl;
        std::cout << "BCtype : " << BCTypes[0] << " , " << BCTypes[1]
                  << std::endl;
        break;
      case 2:
        std::cout << "\\noindent Domain: "
                  << "$(" << lo[0] << "," << (hi[0] - lo[0]) * dx[0]
                  << ")\\times(" << lo[1] << "," << (hi[1] - lo[1]) * dx[1]
                  << ")"
                  << "$\\\\" << std::endl;
        std::cout << "BCtype : " << BCTypes[0] << " , " << BCTypes[1] << " , "
                  << BCTypes[2] << " , " << BCTypes[3] << std::endl;
        break;
    }
    std::cout << multigridID << std::endl;
    std::cout << "\\begin{table}[htbp]" << std::endl;
    std::cout << "\\centering\\begin{tabular}{c|";
    for (int i = 1; i <= 2 * numEncryp + 1; i++)
      std::cout << "c";
    std::cout << "}" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$n$&" << n;
    for (int i = 1; i <= numEncryp; i++) {
      n *= 2;
      std::cout << "&ratio&" << n;
    }
    std::cout << "\\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$\\|\\mathrm{E}\\|_1$&" << res2[0];
    for (int i = 1; i < 2 * numEncryp + 1; i++) {
      std::cout << "&" << res2[i];
    }
    std::cout << "\\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$\\|\\mathrm{E}\\|_2$&" << res3[0];
    for (int i = 1; i < 2 * numEncryp + 1; i++) {
      std::cout << "&" << res3[i];
    }
    std::cout << "\\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$\\|\\mathrm{E}\\|_{\\infty}$&" << res1[0];
    for (int i = 1; i < 2 * numEncryp + 1; i++) {
      std::cout << "&" << res1[i];
    }
    std::cout << "\\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\caption{}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
  }
}

#endif  // _TESTMULTIGRID_H_
