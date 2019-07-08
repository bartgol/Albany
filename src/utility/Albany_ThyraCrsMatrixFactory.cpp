#include "Albany_ThyraCrsMatrixFactory.hpp"

#ifdef ALBANY_EPETRA
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Export.h"
#endif
#include "Albany_TpetraTypes.hpp"

#include "Albany_TpetraThyraUtils.hpp"
#include "Albany_EpetraThyraUtils.hpp"

#include "Albany_Session.hpp"
#include "Albany_Macros.hpp"

namespace Albany {

// The implementation of the graph
struct ThyraCrsMatrixFactory::Impl {

  Impl () = default;

#ifdef ALBANY_EPETRA
  Teuchos::RCP<Epetra_CrsGraph> e_graph;
#endif
  Teuchos::RCP<Tpetra_CrsGraph> t_graph;
};

ThyraCrsMatrixFactory::
ThyraCrsMatrixFactory (const Teuchos::RCP<const Thyra_VectorSpace> domain_vs,
                    const Teuchos::RCP<const Thyra_VectorSpace> range_vs,
                    const int nonzeros_per_row,
                    const bool static_profile)
 : m_graph(new Impl())
 , m_domain_vs(domain_vs)
 , m_range_vs(range_vs)
 , m_filled (false)
 , m_static_profile (static_profile) 
{
  auto bt = Albany::Session::get_build_type();
  TEUCHOS_TEST_FOR_EXCEPTION (bt==BuildType::Undefined, std::logic_error, "Error! No build type set for albany.\n");

  if (bt==BuildType::Epetra) {
#ifdef ALBANY_EPETRA
    auto e_domain = getEpetraBlockMap(domain_vs);
    auto e_range  = getEpetraBlockMap(range_vs);
    m_graph->e_graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*e_range,nonzeros_per_row,static_profile));
#else
    TEUCHOS_TEST_FOR_EXCEPTION (true, std::logic_error, "Error! Epetra is not enabled in albany.\n");
#endif
  } else {
    auto t_domain = getTpetraMap(domain_vs);
    auto t_range = getTpetraMap(range_vs);
    if (static_profile) {
      m_graph->t_graph = Teuchos::rcp(new Tpetra_CrsGraph(t_range,nonzeros_per_row,Tpetra::StaticProfile));
    } else {
      m_graph->t_graph = Teuchos::rcp(new Tpetra_CrsGraph(t_range,nonzeros_per_row,Tpetra::DynamicProfile));
    }
  }
}

ThyraCrsMatrixFactory::
ThyraCrsMatrixFactory (const Teuchos::RCP<const Thyra_VectorSpace> domain_vs,
                       const Teuchos::RCP<const Thyra_VectorSpace> range_vs,
                       const Teuchos::RCP<const ThyraCrsMatrixFactory> overlap_src)
 : m_domain_vs(domain_vs)
 , m_range_vs(range_vs)
{
  TEUCHOS_TEST_FOR_EXCEPTION (!overlap_src->is_filled(), std::logic_error,
                              "Error! Can only build a graph from an overlapped source if source has been filled already.\n");
  m_graph = Teuchos::rcp(new Impl());

  auto bt = Albany::Session::get_build_type();
  TEUCHOS_TEST_FOR_EXCEPTION (bt==BuildType::Undefined, std::logic_error, "Error! No build type set for albany.\n");
  if (bt==BuildType::Epetra) {
#ifdef ALBANY_EPETRA
    auto e_range = getEpetraBlockMap(range_vs);
    auto e_overlap_range = getEpetraBlockMap(overlap_src->m_range_vs);
    auto e_overlap_graph = overlap_src->m_graph->e_graph;

    m_graph->e_graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*e_range,e_overlap_graph->GlobalMaxNumIndices()));

    Epetra_Export exporter (*e_overlap_range,*e_range);
    m_graph->e_graph->Export(*e_overlap_graph,exporter,Insert);

    auto e_domain = getEpetraBlockMap(domain_vs);
    m_graph->e_graph->FillComplete(*e_domain,*e_range);
    m_graph->e_graph->OptimizeStorage();
#else
    TEUCHOS_TEST_FOR_EXCEPTION (true, std::logic_error, "Error! Epetra is not enabled in albany.\n");
#endif
  } else {
    auto t_range = getTpetraMap(range_vs);
    auto t_overlap_range = getTpetraMap(overlap_src->m_range_vs);
    auto t_overlap_graph = overlap_src->m_graph->t_graph;

    m_graph->t_graph = Teuchos::rcp(new Tpetra_CrsGraph(t_range,t_overlap_graph->getGlobalMaxNumRowEntries()));

    Tpetra_Export exporter(t_overlap_range,t_range);
    m_graph->t_graph->doExport(*t_overlap_graph,exporter,Tpetra::INSERT);

    auto t_domain = getTpetraMap(domain_vs);
    auto t_overlap_domain = getTpetraMap(overlap_src->m_domain_vs);
    m_graph->t_graph->fillComplete(t_domain,t_range);
  }

  m_filled = true;
}

void ThyraCrsMatrixFactory::insertGlobalIndices (const GO row, const Teuchos::ArrayView<const GO>& indices) {
  auto bt = Albany::Session::get_build_type();
  if (bt==BuildType::Epetra) {
#ifdef ALBANY_EPETRA
    const GO max_safe_gid = static_cast<GO>(Teuchos::OrdinalTraits<Epetra_GO>::max());
    ALBANY_EXPECT(row<=max_safe_gid, "Error! Input gids exceed Epetra_GO ranges.\n");
    (void) max_safe_gid;

    // Epetra expects pointers to non-const, and Epetra_GO may differ from GO.
    const Epetra_GO e_row = static_cast<Epetra_GO>(row);
    const int e_size = indices.size();
    if (sizeof(GO)==sizeof(Epetra_GO)) {
      // Same size, potentially different type name. A reinterpret_cast will do.
      Epetra_GO* e_indices = const_cast<Epetra_GO*>(reinterpret_cast<const Epetra_GO*>(indices.getRawPtr()));
      m_graph->e_graph->InsertGlobalIndices(e_row, e_size,e_indices);
    } else {
      // Cannot reinterpret cast. Need to copy gids into Epetra_GO array
      Teuchos::Array<Epetra_GO> e_indices(indices.size());
      for (int i=0; i<indices.size(); ++i) {
        ALBANY_EXPECT(indices[i]<=max_safe_gid, "Error! Input gids exceed Epetra_GO ranges.\n");
        e_indices[i] = static_cast<Epetra_GO>(indices[i]);
      }
      m_graph->e_graph->InsertGlobalIndices(e_row, e_size,e_indices.getRawPtr());
    }
#else
    TEUCHOS_TEST_FOR_EXCEPTION (true, std::logic_error, "Error! Epetra is not enabled in albany.\n");
#endif
  } else {
    // Despite being both 64 bits, GO and Tpetra_GO *may* be different *types*.
    Teuchos::ArrayView<const Tpetra_GO> t_indices(reinterpret_cast<const Tpetra_GO*>(indices.getRawPtr()),indices.size());
    m_graph->t_graph->insertGlobalIndices(static_cast<Tpetra_GO>(row),t_indices);
  }
}

void ThyraCrsMatrixFactory::fillComplete () {
  auto bt = Albany::Session::get_build_type();
  if (bt==BuildType::Epetra) {
#ifdef ALBANY_EPETRA
    auto e_domain = getEpetraBlockMap(m_domain_vs);
    auto e_range  = getEpetraBlockMap(m_range_vs);
    m_graph->e_graph->FillComplete(*e_domain,*e_range);
    m_graph->e_graph->OptimizeStorage();
#else
    TEUCHOS_TEST_FOR_EXCEPTION (true, std::logic_error, "Error! Epetra is not enabled in albany.\n");
#endif
  } else {
    auto t_domain = getTpetraMap(m_domain_vs);
    auto t_range = getTpetraMap(m_range_vs);
    m_graph->t_graph->fillComplete(t_domain,t_range);
  }

  m_filled = true;
}

Teuchos::RCP<Thyra_LinearOp> ThyraCrsMatrixFactory::createOp () const {
  TEUCHOS_TEST_FOR_EXCEPTION (!m_filled, std::logic_error, "Error! Cannot create a linear operator if the graph is not filled.\n");

  Teuchos::RCP<Thyra_LinearOp> op;
  auto bt = Albany::Session::get_build_type();
  if (bt==BuildType::Epetra) {
#ifdef ALBANY_EPETRA
    Teuchos::RCP<Epetra_CrsMatrix> matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *m_graph->e_graph)); 
    matrix->PutScalar(0.0); 
    Teuchos::RCP<Epetra_Operator> mat = Teuchos::rcp_implicit_cast<Epetra_Operator>(matrix); 
    op = createThyraLinearOp(mat);
#else
    TEUCHOS_TEST_FOR_EXCEPTION (true, std::logic_error, "Error! Epetra is not enabled in albany.\n");
#endif
  } else {
    Teuchos::RCP<Tpetra_CrsMatrix> mat = Teuchos::rcp (new Tpetra_CrsMatrix(m_graph->t_graph));
    auto const zero = Teuchos::ScalarTraits<ST>::zero();
    mat->resumeFill(); 
    mat->setAllToScalar(zero);
    mat->fillComplete(); 
    op = createThyraLinearOp(Teuchos::rcp_implicit_cast<Tpetra_Operator>(mat));
  }

  return op;
}

} // namespace Albany
