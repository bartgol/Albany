#ifndef ALBANY_THYRA_UTILS_HPP
#define ALBANY_THYRA_UTILS_HPP

// Get basic Thyra types
#include "Albany_ThyraTypes.hpp"

// Get Teuchos_Comm type
#include "Albany_CommTypes.hpp"

// Get Kokkos types (for the 1d device view)
#include "Albany_KokkosTypes.hpp"

namespace Albany
{

// ========= Vector Spaces utilities ========= //

Teuchos::RCP<const Thyra_VectorSpace>
createLocallyReplicatedVectorSpace (const int size, const Teuchos::RCP<const Teuchos_Comm> comm);

GO getGlobalElement (const Teuchos::RCP<const Thyra_VectorSpace>& vs, const LO lid);
LO getLocalElement  (const Teuchos::RCP<const Thyra_VectorSpace>& vs, const GO gid);

// Remove a set of local components from a vector space
Teuchos::RCP<const Thyra_VectorSpace>
removeComponents (const Teuchos::RCP<const Thyra_VectorSpace>& vs,
                  const Teuchos::ArrayView<const LO>& local_components);
// The complement of the above: the specified components are the ones to keep
Teuchos::RCP<const Thyra_VectorSpace>
createSubspace (const Teuchos::RCP<const Thyra_VectorSpace>& vs,
                const Teuchos::ArrayView<const LO>& subspace_components);

// Create a vector space, given the ids of the space components
Teuchos::RCP<const Thyra_VectorSpace>
createVectorSpace (const Teuchos::RCP<const Teuchos_Comm>& comm,
                   const Teuchos::ArrayView<const GO>& gids);

// ========= Thyra_LinearOp utilities ========= //

// These routines help to manipulate the a Thyra_LinearOp. They are needed
// so we can abstract from the concrete linear algebra package, and rely
// only on the Thyra interfaces.

// Fill related helpers
bool isFillActive (const Teuchos::RCP<const Thyra_LinearOp>& lop);
void resumeFill (const Teuchos::RCP<Thyra_LinearOp>& lop);
void fillComplete (const Teuchos::RCP<Thyra_LinearOp>& lop);

// Entries manipulation helpers
void assign (const Teuchos::RCP<Thyra_LinearOp>& lop, const ST value);
void getDiagonalCopy (const Teuchos::RCP<const Thyra_LinearOp>& lop,
                      Teuchos::RCP<Thyra_Vector>& diag);
void getLocalRowValues (const Teuchos::RCP<Thyra_LinearOp>& lop,
                        const LO lrow,
                        Teuchos::Array<LO>& indices,
                        Teuchos::Array<ST>& values);
void setLocalRowValues (const Teuchos::RCP<Thyra_LinearOp>& lop,
                        const LO lrow,
                        const Teuchos::ArrayView<const LO> indices,
                        const Teuchos::ArrayView<const ST> values);
void addToLocalRowValues (const Teuchos::RCP<Thyra_LinearOp>& lop,
                          const LO lrow,
                          const Teuchos::ArrayView<const LO> indices,
                          const Teuchos::ArrayView<const ST> values);

// Math properties helpers
double computeConditionNumber (const Teuchos::RCP<const Thyra_LinearOp>& lop);

// Get a kokkos compatible object to view the content of the linear op on device
DeviceLocalMatrix<const ST> getDeviceData (const Teuchos::RCP<const Thyra_LinearOp>& lop);
DeviceLocalMatrix<ST>       getNonconstDeviceData (const Teuchos::RCP<Thyra_LinearOp>& lop);

// ========= Thyra_(Multi)Vector utilities ========= //

// Const and nonconst version of a getter of local data in a Thyra vector/multivector
Teuchos::ArrayRCP<ST> getNonconstLocalData (const Teuchos::RCP<Thyra_Vector>& v);
Teuchos::ArrayRCP<const ST> getLocalData (const Teuchos::RCP<const Thyra_Vector>& v);
Teuchos::ArrayRCP<Teuchos::ArrayRCP<ST>> getNonconstLocalData (const Teuchos::RCP<Thyra_MultiVector>& mv);
Teuchos::ArrayRCP<Teuchos::ArrayRCP<const ST>> getLocalData (const Teuchos::RCP<const Thyra_MultiVector>& mv);

DeviceView1d<const ST> getDeviceData (const Teuchos::RCP<const Thyra_Vector>& v);
DeviceView1d<ST>       getNonconstDeviceData (const Teuchos::RCP<Thyra_Vector>& v);

// This is just a utility routine, that mildly extend the update method of Thyra_Vector,
// but does not have the complex signature of the linear_combination method of Thyra_Vector.
// In fact, the update method only allows to do y = y + alpha*x, while often one wants
// to do y = beta*y + alpha*x. The linear_combination method offers that capability,
// but the signature is more cumbersome. Thyra offers also a free function with a lighter
// signature for a linear_combination of 2 vectors, but its name is V_StVpStV, which
// can be a bit arcane. Here we simply wrap that function in one with a nicer name.
// Performs y = y_coeff*y + x_coeff*x;
void scale_and_update (const Teuchos::RCP<Thyra_Vector> y, const ST y_coeff,
                       const Teuchos::RCP<const Thyra_Vector> x, const ST x_coeff);

// Thyra does not offer a 'mean' method in its (multi)vector interface.
// The method 'sum' in Thyra_VectorStdOps already does the sum,
// so here we simply scale by the vector (global) length.
ST mean (const Teuchos::RCP<const Thyra_Vector>& v);
Teuchos::Array<ST> means (const Teuchos::RCP<const Thyra_MultiVector>& mv);

// ======== I/O utilities ========= //

template<typename ThyraObjectType>
void describe (const Teuchos::RCP<const ThyraObjectType>& obj,
               Teuchos::FancyOStream& out,
               const Teuchos::EVerbosityLevel verbLevel);

// ========= Thyra_SpmdXYZ utilities ========== //

// These routines help to manipulate Thyra pointers, casting them to
// proper derived classes, and checking that the result is nonnull.

Teuchos::RCP<const Thyra_SpmdVectorSpace>
getSpmdVectorSpace (const Teuchos::RCP<const Thyra_VectorSpace> vs,
                    const bool throw_on_failure = true);

// ========= Thyra_ProductXYZ utilities ========== //

// These routines help to manipulate Thyra pointers, casting them to
// proper derived classes, and checking that the result is nonnull.

Teuchos::RCP<const Thyra_ProductVectorSpace>
getProductVectorSpace (const Teuchos::RCP<const Thyra_VectorSpace> vs,
                       const bool throw_on_failure = true);

Teuchos::RCP<Thyra_ProductVector>
getProductVector (const Teuchos::RCP<Thyra_Vector> v,
                  const bool throw_on_failure = true);

Teuchos::RCP<const Thyra_ProductVector>
getConstProductVector (const Teuchos::RCP<const Thyra_Vector> v,
                       const bool throw_on_failure = true);

Teuchos::RCP<Thyra_ProductMultiVector>
getProductMultiVector (const Teuchos::RCP<Thyra_MultiVector> mv,
                       const bool throw_on_failure = true);

Teuchos::RCP<const Thyra_ProductMultiVector>
getConstProductMultiVector (const Teuchos::RCP<const Thyra_MultiVector> mv,
                            const bool throw_on_failure = true);

} // namespace Albany

#endif // ALBANY_THYRA_UTILS_HPP
