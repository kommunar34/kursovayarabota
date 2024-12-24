Файл linalg.cpp
#include "MatrixWrapper.hpp"
#include "StepIterator.hpp"
#include "simodo/variable/Module_interface.h"
#include "simodo/variable/VariableSetWrapper.h"
#include "simodo/inout/convert/functions.h"
#include "simodo/inout/format/fmt.h"
#include "simodo/bormental/DrBormental.h"
#ifdef CROSS_WIN
// MinGW related workaround
#define BOOST_DLL_FORCE_ALIAS_INSTANTIATION
#endif
#include <boost/dll/alias.hpp>
#include <functional>
#include <limits>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <cmath>
#include <limits>
using namespace simodo;
using namespace simodo::variable;
using namespace simodo::inout;
namespace
{
 Value _vec_dot(Module_interface * , const VariableSetWrapper & args);
 Value _vec_cross(Module_interface * , const VariableSetWrapper & args);
 Value _vector_projection(Module_interface * , const VariableSetWrapper & 
args);
 Value _vector_cos_angle(Module_interface * , const VariableSetWrapper 
&args);
 Value _vector_mixed_product(Module_interface * , const VariableSetWrapper 
&args);
 Value _matr_mul(Module_interface * , const VariableSetWrapper & args);
 Value _matr_apply(Module_interface * , const VariableSetWrapper & args);
 Value _matr_T(Module_interface * , const VariableSetWrapper & args);
 Value _matr_inv(Module_interface * , const VariableSetWrapper & args);
 Value _matr_det(Module_interface * , const VariableSetWrapper & args);
 Value _matr_adj(Module_interface * , const VariableSetWrapper & args);
 Value _matr_id(Module_interface * , const VariableSetWrapper & args);
 
 Value _zero_matrix(Module_interface * , const VariableSetWrapper & args);
 Value _identity_matrix(Module_interface * , const VariableSetWrapper & args);
 Value _matr_pow(Module_interface * , const VariableSetWrapper &args);
 Value _matr_rank(Module_interface * , const VariableSetWrapper &args);
 Value _solve_linear_system(Module_interface * , const VariableSetWrapper 
&args);
 Value _gram_matrix(Module_interface * , const VariableSetWrapper &args);
 Value _vector_norm(Module_interface * , const VariableSetWrapper &args);
 Value _eigenvalues(Module_interface * , const VariableSetWrapper &args);
 Value _pinv(Module_interface * , const VariableSetWrapper &args);
 Value _quat_mul(Module_interface * , const VariableSetWrapper & args);
 Value _quat_apply_and_trunc(Module_interface * , const VariableSetWrapper & 
args);
 Value _quat_dot(Module_interface * , const VariableSetWrapper & args);
 Value _quat_cross(Module_interface * , const VariableSetWrapper & args);
 Value _quat_inv(Module_interface * , const VariableSetWrapper & args);
 Value _quat_conj(Module_interface * , const VariableSetWrapper & args);
 Value _quat_img(Module_interface * , const VariableSetWrapper & args);
 Value _matr_from_quat(Module_interface * , const VariableSetWrapper & args);
 Value _quat_from_matr(Module_interface * , const VariableSetWrapper & args);
}
/*!
 * \brief Получение ссылки на глобальное пространство имён 
математических операций и констант
 */
class ModuleHost_math : public Module_interface
{
public:
 virtual version_t version() const override { return lib_version(); }
 virtual Value instantiate(std::shared_ptr<Module_interface> module_host) 
override {
 return {{
 // Векторы
 {u"vec_dot", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _vec_dot}},
 {{}, ValueType::Float},
 {u"_left_vector_", ValueType::Array},
 {u"_right_vector_", ValueType::Array},
 }}}},
 {u"vec_cross", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _vec_cross}},
 {{}, ValueType::Array},
 {u"_left_vector_", ValueType::Array},
 {u"_right_vector_", ValueType::Array},
 }}}},
 {u"vector_projection", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _vector_projection}},
 {{}, ValueType::Float},
 {u"_a_", ValueType::Array},
 {u"_b_", ValueType::Array},
}}}},
{u"vector_cos_angle", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _vector_cos_angle}},
 {{}, ValueType::Float},
 {u"_a_", ValueType::Array},
 {u"_b_", ValueType::Array},
}}}},
{u"vector_mixed_product", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, 
_vector_mixed_product}},
 {{}, ValueType::Float},
 {u"_a_", ValueType::Array},
 {u"_b_", ValueType::Array},
 {u"_c_", ValueType::Array},
}}}},
 // Матрицы
 // {u"matr_from_vec", {ValueType::Function, Object {{
 // {u"@", ExternalFunction {module_host, _matr_from_vec}},
 // {{}, ValueType::Array},
 // {u"_vector_", ValueType::Array},
 // {u"_is_col_", ValueType::Bool},
 // }}}},
 // {u"matr_row", {ValueType::Function, Object {{
 // {u"@", ExternalFunction {module_host, _matr_row}},
 // {{}, ValueType::Array},
 // {u"_matrix_", ValueType::Array},
 // {u"_row_no_", ValueType::Int},
 // }}}},
 // {u"matr_col", {ValueType::Function, Object {{
 // {u"@", ExternalFunction {module_host, _matr_col}},
 // {{}, ValueType::Array},
 // {u"_matrix_", ValueType::Array},
 // {u"_col_no_", ValueType::Int},
 // }}}},
 {u"pinv", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _pinv}},
 {{}, ValueType::Array},
 {u"_matrix_", ValueType::Array},
}}}},
 {u"eigenvalues", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _eigenvalues}},
 {{}, ValueType::Array},
 {u"_matrix_", ValueType::Array},
}}}},
 {u"vector_norm", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _vector_norm}},
 {{}, ValueType::Float},
 {u"_input_vector_", ValueType::Array},
 {u"_norm_type_", ValueType::Int},
}}}},
 {u"gram_matrix", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _gram_matrix}},
 {{}, ValueType::Array},
 {u"_input_matrix_", ValueType::Array},
}}}},
 {u"solve_linear_system", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, 
_solve_linear_system}},
 {{}, ValueType::Array},
 {u"_matrix_", ValueType::Array},
 {u"_rhs_vector_", ValueType::Array},
}}}},
 {u"matr_rank", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _matr_rank}},
 {{}, ValueType::Int},
 {u"_matrix_", ValueType::Array},
 }}}},
 {u"matr_pow", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _matr_pow}},
 {{}, ValueType::Array},
 {u"_matrix_", ValueType::Array},
 {u"_exponent_", ValueType::Int},
}}}},
 {u"zero_matrix", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _zero_matrix}},
 {{}, ValueType::Array},
 {u"_rows_", ValueType::Int},
 {u"_columns_", ValueType::Int},
 }}}},
 {u"identity_matrix", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _identity_matrix}},
 {{}, ValueType::Array},
 {u"_size_", ValueType::Int},
 }}}},
 {u"matr_mul", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _matr_mul}},
 {{}, ValueType::Array},
 {u"_left_matrix_", ValueType::Array},
 {u"_right_matrix_", ValueType::Array},
 }}}},
 // matr_mul . matr_T $ matr_a, matr_T(vec_b)
 {u"matr_apply", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _matr_apply}},
 {{}, ValueType::Array},
 {u"_matrix_", ValueType::Array},
 {u"_vector_", ValueType::Array},
 }}}},
 {u"matr_T", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _matr_T}},
 {{}, ValueType::Array},
 {u"_matrix_", ValueType::Array},
 }}}},
 {u"matr_inv", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _matr_inv}},
 {{}, ValueType::Array},
 {u"_matrix_", ValueType::Array},
 }}}},
 {u"matr_det", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _matr_det}},
 {{}, ValueType::Float},
 {u"_matrix_", ValueType::Array},
 }}}},
 {u"matr_adj", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _matr_adj}},
 {{}, ValueType::Array},
 {u"_matrix_", ValueType::Array},
 }}}},
 {u"matr_id", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _matr_id}},
 {{}, ValueType::Array},
 {u"_side_size_", ValueType::Int},
 }}}},
 // Кватернионы
 {u"quat_mul", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _quat_mul}},
 {{}, ValueType::Array},
 {u"_left_quaternion_", ValueType::Array},
 {u"_right_quaternion_", ValueType::Array},
 }}}},
 {u"quat_apply_and_trunc", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _quat_apply_and_trunc}},
 {{}, ValueType::Array},
 {u"_quaternion_", ValueType::Array},
 {u"_vector_", ValueType::Array},
 }}}},
 {u"quat_dot", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _quat_dot}},
 {{}, ValueType::Float},
 {u"_left_quaternion_", ValueType::Array},
 {u"_right_quaternion_", ValueType::Array},
 }}}},
 {u"quat_cross", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _quat_cross}},
 {{}, ValueType::Array},
 {u"_left_quaternion_", ValueType::Array},
 {u"_right_quaternion_", ValueType::Array},
 }}}},
 {u"quat_inv", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _quat_inv}},
 {{}, ValueType::Array},
 {u"_quaternion_", ValueType::Array},
 }}}},
 {u"quat_conj", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _quat_conj}},
 {{}, ValueType::Array},
 {u"_quaternion_", ValueType::Array},
 }}}},
 {u"quat_img", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _quat_img}},
 {{}, ValueType::Array},
 {u"_quaternion_", ValueType::Array},
 }}}},
 {u"matr_from_quat", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _matr_from_quat}},
 {{}, ValueType::Array},
 {u"_quaternion_", ValueType::Array},
 }}}},
 {u"quat_from_matr", {ValueType::Function, Object {{
 {u"@", ExternalFunction {module_host, _quat_from_matr}},
 {{}, ValueType::Array},
 {u"_rotation_matrix_", ValueType::Array},
 }}}},
 }};
 }
 // Factory method
 static std::shared_ptr<Module_interface> create() {
 return std::make_shared<ModuleHost_math>();
 }
};
BOOST_DLL_ALIAS(
 ModuleHost_math::create, // <-- this function is exported with...
 create_simodo_module // <-- ...this alias name
)
namespace
{
 bool is_valid_array(const Array &array)
 {
 const auto &dimensions = array.dimensions();
 constexpr auto index_is_positive
 = [](index_t i) -> bool { return i > 0; };
 return ((dimensions.empty()
 || (dimensions.size() == 1
 && dimensions[0] == 0))
 && array.values().empty())
 || (std::all_of(
 dimensions.cbegin(),
 dimensions.cend(),
 index_is_positive
 )
 && std::reduce(
 dimensions.cbegin(),
 dimensions.cend(),
 uint16_t(1),
 std::multiplies<>()
 ) == array.values().size());
 }
 void assert_number(const Value &value)
 {
 if (value.isNumber()) return;
 throw bormental::DrBormental(
 "Module linalg assert_number", 
 inout::fmt("Value with type %1 is not a number.")
 .arg(getValueTypeName(value.type())));
 }
 double get_number(const Value &value)
 {
 if (value.isInt()) return value.getInt();
 return value.getFloat();
 }
 Value plus_values(const Value &lhs, const Value &rhs)
 {
 assert_number(lhs);
 assert_number(rhs);
 if (lhs.isInt() && rhs.isInt())
 {
 return int64_t(lhs.getInt() + rhs.getInt());
 }
 return double(get_number(lhs) + get_number(rhs));
 }
 Value minus_values(const Value &lhs, const Value &rhs)
 {
 assert_number(lhs);
 assert_number(rhs);
 if (lhs.isInt() && rhs.isInt())
 {
 return int64_t(lhs.getInt() - rhs.getInt());
 }
 return double(get_number(lhs) - get_number(rhs));
 }
 Value multiply_values(const Value &lhs, const Value &rhs)
 {
 assert_number(lhs);
 assert_number(rhs);
 if (lhs.isInt() && rhs.isInt())
 {
 return int64_t(lhs.getInt() * rhs.getInt());
 }
 return double(get_number(lhs) * get_number(rhs));
 }
 Value _vec_dot(Module_interface * , const VariableSetWrapper & args)
 {
 const Value & lhs = args[0].origin().value();
 assert(lhs.isArray() && is_valid_array(*lhs.getArray()));
 const Value & rhs = args[1].origin().value();
 assert(rhs.isArray() && is_valid_array(*rhs.getArray()));
 const auto lhs_array = lhs.getArray();
 const auto rhs_array = rhs.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 const auto &rhs_dimensions = rhs_array->dimensions();
 if (lhs_dimensions.size() != 1 || lhs_dimensions.size() != 
rhs_dimensions.size())
 {
 throw bormental::DrBormental(
 "Module linalg vec_dot",
 "The arguments are not both the same length 1D array[*].");
 }
 if (lhs_dimensions[0] != rhs_dimensions[0])
 {
 throw bormental::DrBormental(
 "Module linalg vec_dot",
 fmt("The arguments are not both"
 " same length 1D array[*]: left = %1,"
 " right = %2.")
 .arg(lhs_dimensions[0])
 .arg(rhs_dimensions[0])
 );
 }
 return std::transform_reduce(
 lhs_array->values().cbegin(),
 lhs_array->values().cend(),
 rhs_array->values().cbegin(),
 Value(int64_t(0)),
 plus_values,
 multiply_values
 );
 }
 Value _vec_cross(Module_interface * , const VariableSetWrapper & args)
 {
 const Value & lhs_value = args[0].origin().value();
 assert(lhs_value.isArray() && is_valid_array(*lhs_value.getArray()));
 const Value & rhs_value = args[1].origin().value();
 assert(rhs_value.isArray() && is_valid_array(*rhs_value.getArray()));
 const auto &lhs_array = lhs_value.getArray();
 const auto &rhs_array = rhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 const auto &rhs_dimensions = rhs_array->dimensions();
 if (lhs_dimensions.size() != 1 || lhs_dimensions.size() != 
rhs_dimensions.size())
 {
 throw bormental::DrBormental(
 "Module linalg vec_cross",
 "The arguments are not both 1D array[3].");
 }
 if (lhs_dimensions[0] != 3 || lhs_dimensions[0] != rhs_dimensions[0])
 {
 throw bormental::DrBormental(
 "Module linalg vec_cross",
 fmt("The arguments are not both 1D arrays[3]: left = %1,"
 " right = %2.")
 .arg(lhs_dimensions[0])
 .arg(rhs_dimensions[0])
 );
 }
 const auto &lhs = lhs_array->values();
 const auto &rhs = rhs_array->values();
 return std::make_shared<Array>(std::vector<Value>{
 minus_values(
 multiply_values(lhs[1], rhs[2]),
 multiply_values(lhs[2], rhs[1])
 ),
 minus_values(
 multiply_values(lhs[2], rhs[0]),
 multiply_values(lhs[0], rhs[2])
 ),
 minus_values(
 multiply_values(lhs[0], rhs[1]),
 multiply_values(lhs[1], rhs[0])
 ),
 });
 }
 Value _matr_mul(Module_interface * , const VariableSetWrapper & args)
 {
 const Value & lhs_value = args[0].origin().value();
 assert(lhs_value.isArray() && is_valid_array(*lhs_value.getArray()));
 const Value & rhs_value = args[1].origin().value();
 assert(rhs_value.isArray() && is_valid_array(*rhs_value.getArray()));
 const auto &lhs_array = lhs_value.getArray();
 const auto &rhs_array = rhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 const auto &rhs_dimensions = rhs_array->dimensions();
 if (lhs_dimensions.size() != 2 || lhs_dimensions.size() != 
rhs_dimensions.size())
 {
 throw bormental::DrBormental(
 "Module linalg matr_mul",
 "The arguments are not both 2D array[n,k] and array[k,m].");
 }
 if (lhs_dimensions[1] != rhs_dimensions[0])
 {
 throw bormental::DrBormental(
 "Module linalg matr_mul",
 fmt("The arguments are not both 2D array[n,k] and array[k,m]:"
 " left = [*,%1], right = [%2,*].")
 .arg(lhs_dimensions[1])
 .arg(rhs_dimensions[0])
 );
 }
 const auto &lhs = lhs_array->values();
 const auto &rhs = rhs_array->values();
 index_t n_product_rows = lhs_dimensions[0];
 index_t n_product_cols = rhs_dimensions[1];
 auto product = std::make_shared<Array>(
 ValueType::Null,
 std::vector<index_t>{
 index_t(n_product_rows),
 index_t(n_product_cols),
 },
 std::vector<Value>()
 );
 product->values_mutable().reserve(n_product_rows * n_product_cols);
 for (index_t i = 0; i < n_product_rows; ++i)
 {
 for (index_t j = 0; j < n_product_cols; ++j)
 {
 auto row_begin
 = lhs.cbegin() + lhs_dimensions[1] * i;
 auto row_end
 = lhs.cbegin() + lhs_dimensions[1] * (i + 1);
 auto col_begin = StepIterator< Value >(
 &rhs[j],
 typename StepIterator<
 Value >::difference_type(n_product_cols)
 );
 auto init_value = Value(int64_t(0));
 product->add(std::transform_reduce(
 row_begin,
 row_end,
 col_begin,
 init_value,
 plus_values,
 multiply_values
 ));
 }
 }
 return product;
 }
 Value _matr_apply(
 Module_interface *, const VariableSetWrapper &args
 )
 {
 const Value &lhs_value = args[0].origin().value();
 assert(
 lhs_value.isArray()
 && is_valid_array(*lhs_value.getArray())
 );
 const Value &rhs_value = args[1].origin().value();
 assert(
 rhs_value.isArray()
 && is_valid_array(*rhs_value.getArray())
 );
 const auto &lhs_array = lhs_value.getArray();
 const auto &rhs_array = rhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 const auto &rhs_dimensions = rhs_array->dimensions();
 if (lhs_dimensions.size() != 2
 || rhs_dimensions.size() != 1)
 {
 throw bormental::DrBormental(
 "Module linalg matr_apply",
 "The left argument is not 2D array[n,k] and "
 "the right argument is not 1D array[k]."
 );
 }
 if (lhs_dimensions[1] != rhs_dimensions[0])
 {
 throw bormental::DrBormental(
 "Module linalg matr_apply",
 fmt("The left argument is not 2D array[n,k] "
 "or the right argument is not 1D array[k]:"
 " left = [*,%1], right = [%2].")
 .arg(lhs_dimensions[1])
 .arg(rhs_dimensions[0])
 );
 }
 const auto &lhs = lhs_array->values();
 const auto &rhs = rhs_array->values();
 index_t n_product_rows = lhs_dimensions[0];
 std::vector< Value > product;
 product.reserve(n_product_rows);
 for (index_t i = 0; i < n_product_rows; ++i)
 {
 auto row_begin
 = lhs.cbegin() + lhs_dimensions[1] * i;
 auto row_end
 = lhs.cbegin() + lhs_dimensions[1] * (i + 1);
 auto col_begin = rhs.cbegin();
 auto init_value = Value(int64_t(0));
 product.push_back(std::transform_reduce(
 row_begin,
 row_end,
 col_begin,
 init_value,
 plus_values,
 multiply_values
 ));
 }
 return std::make_shared< Array >(std::move(product));
 }
 Value _matr_T(
 Module_interface *, const VariableSetWrapper &args
 )
 {
 const Value &lhs_value = args[0].origin().value();
 assert(
 lhs_value.isArray()
 && is_valid_array(*lhs_value.getArray())
 );
 const auto &lhs_array = lhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 if (lhs_dimensions.size() != 2)
 {
 throw bormental::DrBormental(
 "Module linalg matr_T",
 "The argument is not 2D array[*,*]."
 );
 }
 const auto &lhs = lhs_array->values();
 index_t n_transposed_rows = lhs_dimensions[0];
 index_t n_transposed_cols = lhs_dimensions[1];
 std::vector< Value > transposed;
 transposed.reserve(
 n_transposed_rows * n_transposed_cols
 );
 for (index_t i = 0; i < n_transposed_rows; ++i)
 {
 auto col_begin = StepIterator< Value >(
 &lhs[i],
 typename StepIterator< Value >::difference_type(
 n_transposed_rows
 )
 );
 auto col_end = col_begin + n_transposed_cols;
 std::copy(
 col_begin,
 col_end,
 std::back_inserter(transposed)
 );
 }
 return std::make_shared< Array >(
 ValueType::Null,
 std::vector< index_t >{
 index_t(n_transposed_rows),
 index_t(n_transposed_cols),
 },
 std::move(transposed)
 );
 }
 template < typename T > inline bool real_is_zero(T t)
 {
 return std::abs(t)
 <= std::abs(t) * 2
 * std::numeric_limits< T >::epsilon();
 }
 Value _matr_inv(
 Module_interface *, const VariableSetWrapper &args
 )
 {
 const Value &lhs_value = args[0].origin().value();
 assert(
 lhs_value.isArray()
 && is_valid_array(*lhs_value.getArray())
 );
 const auto &lhs_array = lhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 if (lhs_dimensions.size() != 2)
 {
 throw bormental::DrBormental(
 "Module linalg matr_inv",
 "The argument is not 2D array[k,k]."
 );
 }
 if (lhs_dimensions[0] != lhs_dimensions[1])
 {
 throw bormental::DrBormental(
 "Module linalg matr_inv",
 fmt("The argument is not 2D array[k,k]: "
 "[%1,%2]")
 .arg(lhs_dimensions[0])
 .arg(lhs_dimensions[1])
 );
 }
 index_t matrix_side = lhs_dimensions[0];
 MatrixWrapper< Value, index_t > matrix_wrapper
 = {matrix_side, lhs_array->values()};
 Value determinant = matrix_wrapper.determinant(
 plus_values, minus_values, multiply_values
 );
 bool is_matrix_invertible
 = (determinant.isInt())
 ? (determinant.getInt() != 0)
 : (!real_is_zero< double >(determinant.getFloat(
 )));
 if (!is_matrix_invertible)
 {
 throw bormental::DrBormental(
 "Module linalg matr_inv",
 "The argument is not invertible matrix."
 );
 }
 Value inversed_determinant = double(
 1.
 / (determinant.isInt()
 ? double(determinant.getInt())
 : determinant.getFloat())
 );
 std::vector< Value > inversed
 = matrix_wrapper.transposed_cofactor_matrix(
 plus_values, minus_values, multiply_values
 );
 std::transform(
 inversed.cbegin(),
 inversed.cend(),
 inversed.begin(),
 [inversed_determinant](const Value &v)
 { return multiply_values(inversed_determinant, v); }
 );
 return std::make_shared< Array >(
 ValueType::Null,
 std::vector< index_t >{
 index_t(matrix_side),
 index_t(matrix_side),
 },
 std::move(inversed)
 );
 }
 Value _matr_det(Module_interface * , const VariableSetWrapper & args)
 {
 const Value &lhs_value = args[0].origin().value();
 assert(
 lhs_value.isArray()
 && is_valid_array(*lhs_value.getArray())
 );
 const auto &lhs_array = lhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 if (lhs_dimensions.size() != 2)
 {
 throw bormental::DrBormental(
 "Module linalg matr_det",
 "The argument is not 2D array[k,k]."
 );
 }
 if (lhs_dimensions[0] != lhs_dimensions[1])
 {
 throw bormental::DrBormental(
 "Module linalg matr_det",
 fmt("The argument is not 2D array[k,k]: "
 "[%1,%2]")
 .arg(lhs_dimensions[0])
 .arg(lhs_dimensions[1])
 );
 }
 index_t matrix_side = lhs_dimensions[0];
 MatrixWrapper< Value, index_t > matrix_wrapper
 = {matrix_side, lhs_array->values()};
 return matrix_wrapper.determinant(
 plus_values, minus_values, multiply_values
 );
 }
 Value _matr_adj(Module_interface * , const VariableSetWrapper & args)
 {
 const Value &lhs_value = args[0].origin().value();
 assert(
 lhs_value.isArray()
 && is_valid_array(*lhs_value.getArray())
 );
 const auto &lhs_array = lhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 if (lhs_dimensions.size() != 2)
 {
 throw bormental::DrBormental(
 "Module linalg matr_adj",
 "The argument is not 2D array[k,k]."
 );
 }
 if (lhs_dimensions[0] != lhs_dimensions[1])
 {
 throw bormental::DrBormental(
 "Module linalg matr_adj",
 fmt("The argument is not 2D array[k,k]: "
 "[%1,%2]")
 .arg(lhs_dimensions[0])
 .arg(lhs_dimensions[1])
 );
 }
 index_t matrix_side = lhs_dimensions[0];
 MatrixWrapper< Value, index_t > matrix_wrapper
 = {matrix_side, lhs_array->values()};
 return std::make_shared< Array >(
 ValueType::Null,
 std::vector< index_t >{
 index_t(matrix_side),
 index_t(matrix_side),
 },
 matrix_wrapper.transposed_cofactor_matrix(
 plus_values, minus_values, multiply_values
 )
 );
 }
 Value _matr_id(
 Module_interface *, const VariableSetWrapper &args
 )
 {
 const Value &lhs_value = args[0].origin().value();
 assert(lhs_value.isInt());
 int64_t lhs = lhs_value.getInt();
 if (lhs < 1 || lhs > 32'654)
 {
 throw bormental::DrBormental(
 "Module linalg matr_id",
 fmt("The argument must in [1;32'654]: %1")
 .arg(lhs)
 );
 }
 std::vector< Value > identity;
 identity.reserve(lhs * lhs);
 for (index_t i = 0; i < lhs; ++i)
 {
 for (index_t j = 0; j < lhs; ++j)
 {
 identity.emplace_back(
 i == j ? int64_t(1) : int64_t(0)
 );
 }
 }
 return std::make_shared< Array >(
 ValueType::Null,
 std::vector< index_t >{
 index_t(lhs),
 index_t(lhs),
 },
 std::move(identity)
 );
 }
 template <
 typename T,
 typename I,
 typename C = std::vector< T > >
 struct hamilton_product
 {
 template <
 typename C1,
 typename C2,
 typename P,
 typename Mi,
 typename M >
 C operator()(
 const C1 &lhs,
 const C2 &rhs,
 P plus,
 Mi minus,
 M mul
 ) const
 {
 T dot_product = std::transform_reduce(
 lhs.cbegin() + 1,
 lhs.cend(),
 rhs.cbegin() + 1,
 T(I(0)),
 plus,
 mul
 );
 T real_part
 = minus(mul(lhs[0], rhs[0]), dot_product);
 std::array< T, 3 > imaginary_part = {
 minus(mul(lhs[2], rhs[3]), mul(lhs[3], rhs[2])),
 minus(mul(lhs[3], rhs[1]), mul(lhs[1], rhs[3])),
 minus(mul(lhs[1], rhs[2]), mul(lhs[2], rhs[1])),
 };
 std::transform(
 imaginary_part.cbegin(),
 imaginary_part.cend(),
 rhs.cbegin() + 1,
 imaginary_part.begin(),
 [plus,
 mul,
 &lhs](const T &img_v, const T &rhs_v) -> T
 { return plus(img_v, mul(lhs[0], rhs_v)); }
 );
 std::transform(
 imaginary_part.cbegin(),
 imaginary_part.cend(),
 lhs.cbegin() + 1,
 imaginary_part.begin(),
 [plus,
 mul,
 &rhs](const T &img_v, const T &lhs_v) -> T
 { return plus(img_v, mul(rhs[0], lhs_v)); }
 );
 return C{
 real_part,
 imaginary_part[0],
 imaginary_part[1],
 imaginary_part[2],
 };
 }
 
 
 };
 Value _quat_mul(
 Module_interface *, const VariableSetWrapper &args
 )
 {
 const Value &lhs_value = args[0].origin().value();
 assert(
 lhs_value.isArray()
 && is_valid_array(*lhs_value.getArray())
 );
 const Value &rhs_value = args[1].origin().value();
 assert(
 rhs_value.isArray()
 && is_valid_array(*rhs_value.getArray())
 );
 const auto &lhs_array = lhs_value.getArray();
 const auto &rhs_array = rhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 const auto &rhs_dimensions = rhs_array->dimensions();
 if (lhs_dimensions.size() != 1
 || lhs_dimensions.size() != rhs_dimensions.size())
 {
 throw bormental::DrBormental(
 "Module linalg quat_mul",
 "The arguments are not both 1D array[4]."
 );
 }
 if (lhs_dimensions[0] != 4
 || lhs_dimensions[0] != rhs_dimensions[0])
 {
 throw bormental::DrBormental(
 "Module linalg quat_mul",
 fmt("The arguments are not both 1D arrays[4]: "
 "left = %1, right = %2.")
 .arg(lhs_dimensions[0])
 .arg(rhs_dimensions[0])
 );
 }
 const auto &lhs = lhs_array->values();
 const auto &rhs = rhs_array->values();
 return std::make_shared<
 Array >(hamilton_product< Value, int64_t >{}(
 lhs, rhs, plus_values, minus_values, multiply_values
 ));
 }
 Value _quat_apply_and_trunc(
 Module_interface *, const VariableSetWrapper &args
 )
 {
 const Value &lhs_value = args[0].origin().value();
 assert(
 lhs_value.isArray()
 && is_valid_array(*lhs_value.getArray())
 );
 const Value &rhs_value = args[1].origin().value();
 assert(
 rhs_value.isArray()
 && is_valid_array(*rhs_value.getArray())
 );
 const auto &lhs_array = lhs_value.getArray();
 const auto &rhs_array = rhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 const auto &rhs_dimensions = rhs_array->dimensions();
 if (lhs_dimensions.size() != 1
 || lhs_dimensions.size() != rhs_dimensions.size())
 {
 throw bormental::DrBormental(
 "Module linalg quat_apply_and_trunc",
 "The left argument is not 1D array[4] "
 "or the right argument is not 1D array[3]."
 );
 }
 if (lhs_dimensions[0] != 4 || rhs_dimensions[0] != 3)
 {
 throw bormental::DrBormental(
 "Module linalg quat_apply_and_trunc",
 fmt("The left argument is not 1D array[4] "
 "or the right argument is not 1D array[3]: "
 "left = %1, right = %2.")
 .arg(lhs_dimensions[0])
 .arg(rhs_dimensions[0])
 );
 }
 const auto &q = lhs_array->values();
 std::array<Value, 4> v;
 v[0] = int64_t(0);
 std::copy(
 rhs_array->values().cbegin(),
 rhs_array->values().cend(),
 v.begin() + 1
 );
 hamilton_product< Value, int64_t> ham_prod;
 auto untruncated_result = ham_prod(
 ham_prod(
 q, v, plus_values, minus_values, multiply_values
 ),
 std::vector< Value >{
 q[0],
 minus_values(int64_t(0), q[1]),
 minus_values(int64_t(0), q[2]),
 minus_values(int64_t(0), q[3])
 },
 plus_values,
 minus_values,
 multiply_values
 );
 return std::make_shared< Array >(std::vector< Value >(
 untruncated_result.begin() + 1,
 untruncated_result.end()
 ));
 }
 Value _quat_dot(Module_interface * , const VariableSetWrapper & /*args*/)
 {
 throw bormental::DrBormental(
 "Module linalg quat_dot",
 "Current operation is not supported");
 }
 Value _quat_cross(Module_interface * , const VariableSetWrapper & /*args*/)
 {
 throw bormental::DrBormental(
 "Module linalg quat_cross",
 "Current operation is not supported");
 }
 Value _quat_inv(
 Module_interface *, const VariableSetWrapper &args
 )
 {
 const Value &lhs_value = args[0].origin().value();
 assert(
 lhs_value.isArray()
 && is_valid_array(*lhs_value.getArray())
 );
 const auto &lhs_array = lhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 if (lhs_dimensions.size() != 1)
 {
 throw bormental::DrBormental(
 "Module linalg quat_inv",
 "The argument is not 1D array[4]."
 );
 }
 if (lhs_dimensions[0] != 4)
 {
 throw bormental::DrBormental(
 "Module linalg quat_inv",
 fmt("The argument is not 1D arrays[4]: %1.")
 .arg(lhs_dimensions[0])
 );
 }
 const auto &lhs = lhs_array->values();
 std::array< double, 4 > q;
 std::transform(
 lhs.cbegin(),
 lhs.cend(),
 q.begin(),
 get_number
 );
 double norm = std::sqrt(std::transform_reduce(
 q.cbegin(),
 q.cend(),
 double(0.),
 std::plus<>(),
 [](double v) -> double { return v * v; }
 ));
 if (real_is_zero< double >(norm))
 {
 throw bormental::DrBormental(
 "Module linalg quat_inv",
 "The argument is not invertible quaternion"
 );
 }
 std::transform(
 q.cbegin(),
 q.cend(),
 q.begin(),
 [norm](double v) -> double { return v / norm; }
 );
 return std::make_shared< Array >(std::vector< Value >{
 q[0],
 -q[1],
 -q[2],
 -q[3],
 });
 }
 Value _quat_conj(
 Module_interface *, const VariableSetWrapper &args
 )
 {
 const Value &lhs_value = args[0].origin().value();
 assert(
 lhs_value.isArray()
 && is_valid_array(*lhs_value.getArray())
 );
 const auto &lhs_array = lhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 if (lhs_dimensions.size() != 1)
 {
 throw bormental::DrBormental(
 "Module linalg quat_img",
 "The argument is not 1D array[4]."
 );
 }
 if (lhs_dimensions[0] != 4)
 {
 throw bormental::DrBormental(
 "Module linalg quat_img",
 fmt("The argument is not 1D arrays[4]: %1.")
 .arg(lhs_dimensions[0])
 );
 }
 const auto &lhs = lhs_array->values();
 return std::make_shared< Array >(std::vector< Value >{
 lhs[0],
 minus_values(int64_t(0), lhs[1]),
 minus_values(int64_t(0), lhs[2]),
 minus_values(int64_t(0), lhs[3]),
 });
 }
 Value _quat_img(
 Module_interface *, const VariableSetWrapper &args
 )
 {
 const Value &lhs_value = args[0].origin().value();
 assert(
 lhs_value.isArray()
 && is_valid_array(*lhs_value.getArray())
 );
 const auto &lhs_array = lhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 if (lhs_dimensions.size() != 1)
 {
 throw bormental::DrBormental(
 "Module linalg quat_img",
 "The argument is not 1D array[4]."
 );
 }
 if (lhs_dimensions[0] != 4)
 {
 throw bormental::DrBormental(
 "Module linalg quat_img",
 fmt("The argument is not 1D arrays[4]: %1.")
 .arg(lhs_dimensions[0])
 );
 }
 const auto &lhs = lhs_array->values();
 return std::make_shared< Array >(
 std::vector< Value >{lhs[1], lhs[2], lhs[3]}
 );
 }
 constexpr double
 saturate(double val, double min, double max)
 {
 return std::min(max, std::max(min, val));
 }
 constexpr double saturate_abs_1(double val)
 {
 return saturate(val, -1., 1.);
 }
 // 1. - 2. * (q[2] * q[2] + q[3] * q[3]),
 // 0. + 2. * (q[1] * q[2] - q[3] * q[0]),
 // 0. + 2. * (q[1] * q[3] + q[2] * q[0]),
 using const_quaternion_ref_t
 = const std::array< double, 4 > &;
 inline double _matrX__from_quat(const_quaternion_ref_t q)
 {
 return saturate_abs_1(
 q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3]
 );
 }
 inline double _matrXY_from_quat(const_quaternion_ref_t q)
 {
 return saturate_abs_1(
 0. + 2. * (q[1] * q[2] - q[3] * q[0])
 );
 }
 inline double _matrXZ_from_quat(const_quaternion_ref_t q)
 {
 return saturate_abs_1(
 0. + 2. * (q[1] * q[3] + q[2] * q[0])
 );
 }
 // 0. + 2. * (q[1] * q[2] + q[3] * q[0]),
 // 1. - 2. * (q[1] * q[1] + q[3] * q[3]),
 // 0. + 2. * (q[2] * q[3] - q[1] * q[0]),
 inline double _matrYX_from_quat(const_quaternion_ref_t q)
 {
 return saturate_abs_1(
 0. + 2. * (q[1] * q[2] + q[3] * q[0])
 );
 }
 inline double _matrY__from_quat(const_quaternion_ref_t q)
 {
 return saturate_abs_1(
 q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3]
 );
 }
 inline double _matrYZ_from_quat(const_quaternion_ref_t q)
 {
 return saturate_abs_1(
 0. + 2. * (q[2] * q[3] - q[1] * q[0])
 );
 }
 // 0. + 2. * (q[1] * q[3] - q[2] * q[0]),
 // 0. + 2. * (q[2] * q[3] + q[1] * q[0]),
 // 1. - 2. * (q[1] * q[1] + q[2] * q[2]),
 inline double _matrZX_from_quat(const_quaternion_ref_t q)
 {
 return saturate_abs_1(
 0. + 2. * (q[1] * q[3] - q[2] * q[0])
 );
 }
 inline double _matrZY_from_quat(const_quaternion_ref_t q)
 {
 return saturate_abs_1(
 0. + 2. * (q[2] * q[3] + q[1] * q[0])
 );
 }
 inline double _matrZ__from_quat(const_quaternion_ref_t q)
 {
 return saturate_abs_1(
 q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]
 );
 }
 Value _matr_from_quat(
 Module_interface *, const VariableSetWrapper &args
 )
 {
 const Value &lhs_value = args[0].origin().value();
 assert(
 lhs_value.isArray()
 && is_valid_array(*lhs_value.getArray())
 );
 const auto &lhs_array = lhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 if (lhs_dimensions.size() != 1)
 {
 throw bormental::DrBormental(
 "Module linalg matr_from_quat",
 "The argument is not 1D array[4]."
 );
 }
 if (lhs_dimensions[0] != 4)
 {
 throw bormental::DrBormental(
 "Module linalg matr_from_quat",
 fmt("The argument is not 1D arrays[4]: %1.")
 .arg(lhs_dimensions[0])
 );
 }
 const auto &lhs = lhs_array->values();
 std::array< double, 4 > q;
 std::transform(
 lhs.cbegin(),
 lhs.cend(),
 q.begin(),
 get_number
 );
 // Нормирование
 // double norm = std::sqrt(std::transform_reduce(
 // q.cbegin(),
 // q.cend(),
 // double(0.),
 // std::plus<>(),
 // [](double v) -> double { return v * v; }
 // ));
 // std::transform(
 // q.cbegin(),
 // q.cend(),
 // q.begin(),
 // [norm](double v) -> double { return v / norm; }
 // );
 return std::make_shared< Array >(
 ValueType::Null,
 std::vector< index_t >{index_t(3), index_t(3)},
 std::vector< Value >{
 _matrX__from_quat(q),
 _matrXY_from_quat(q),
 _matrXZ_from_quat(q),
 _matrYX_from_quat(q),
 _matrY__from_quat(q),
 _matrYZ_from_quat(q),
 _matrZX_from_quat(q),
 _matrZY_from_quat(q),
 _matrZ__from_quat(q),
 }
 );
 }
 Value _quat_from_matr(Module_interface * , const VariableSetWrapper & 
/*args*/)
 {
 throw bormental::DrBormental(
 "Module linalg quat_from_matr",
 "Current operation is not supported");
 }
 
 Value _zero_matrix(Module_interface * , const VariableSetWrapper & args) {
 const int rows = args[0].origin().value().getInt();
 const int columns = args[1].origin().value().getInt();
 
 std::vector<Value> values(rows * columns, Value(0.0));
 std::vector<index_t> dimensions = {index_t(rows), index_t(columns)};
 
 return std::make_shared<Array>(ValueType::Null, dimensions, values); 
 }
 
 Value _identity_matrix(Module_interface * , const VariableSetWrapper & args) 
{
 const int size = args[0].origin().value().getInt();
 
 std::vector<Value> values(size * size, Value(0.0));
 for (int i = 0; i < size; ++i) {
 values[i * size + i] = Value(1.0);
 }
 
 std::vector<index_t> dimensions = {index_t(size), index_t(size)};
 
 return std::make_shared<Array>(ValueType::Null, dimensions, values);
 } 
 Value _vector_projection(Module_interface * , const VariableSetWrapper & 
args) 
 {
 const Value & lhs_value = args[0].origin().value();
 assert(lhs_value.isArray() && is_valid_array(*lhs_value.getArray()));
 const Value & rhs_value = args[1].origin().value();
 assert(rhs_value.isArray() && is_valid_array(*rhs_value.getArray()));
 const auto &lhs_array = lhs_value.getArray();
 const auto &rhs_array = rhs_value.getArray();
 const auto &lhs_dimensions = lhs_array->dimensions();
 const auto &rhs_dimensions = rhs_array->dimensions();
 if (lhs_dimensions.size() != 1 || lhs_dimensions.size() != 
rhs_dimensions.size())
 {
 throw bormental::DrBormental(
 "Module linalg vec_cross",
 "The arguments are not both 1D array[3].");
 }
 if (lhs_dimensions[0] != 3 || lhs_dimensions[0] != rhs_dimensions[0])
 {
 throw bormental::DrBormental(
 "Module linalg vec_cross",
 fmt("The arguments are not both 1D arrays[3]: left = %1,"
 " right = %2.")
 .arg(lhs_dimensions[0])
 .arg(rhs_dimensions[0])
 );
 }
 const auto &lhs = lhs_array->values();
 const auto &rhs = rhs_array->values();
 
 double dot_product_ab = 0.0;
 for (size_t i = 0; i < lhs.size(); ++i) {
 dot_product_ab += get_number(lhs[i]) * get_number(rhs[i]);
 }
 double dot_product_bb = 0.0;
 for (size_t i = 0; i < rhs.size(); ++i) {
 dot_product_bb += get_number(rhs[i]) * get_number(rhs[i]);
 }
 double projection_scale = dot_product_ab / std::sqrt(dot_product_bb);
return Value(projection_scale);
 }
Value _vector_cos_angle(Module_interface *, const VariableSetWrapper &args) {
 const Value &a_value = args[0].origin().value();
 const Value &b_value = args[1].origin().value();
 assert(a_value.isArray() && b_value.isArray() &&
 is_valid_array(*a_value.getArray()) && 
is_valid_array(*b_value.getArray()));
 const auto &a_array = a_value.getArray();
 const auto &b_array = b_value.getArray();
 const auto &a_values = a_array->values();
 const auto &b_values = b_array->values();
 // Проверяем, что оба вектора имеют одинаковую размерность
 if (a_values.size() != b_values.size()) {
 throw bormental::DrBormental(
 "Module linalg vector_cos_angle",
 "The vectors must have the same size"
 );
 }
 // Скалярное произведение a и b
 double dot_product_ab = 0.0;
 for (size_t i = 0; i < a_values.size(); ++i) {
 dot_product_ab += get_number(multiply_values(a_values[i], b_values[i]));
 }
 // Норма вектора a
 double norm_a = 0.0;
 for (size_t i = 0; i < a_values.size(); ++i) {
 norm_a += get_number(multiply_values(a_values[i], a_values[i]));
 }
 norm_a = std::sqrt(norm_a);
 // Норма вектора b
 double norm_b = 0.0;
 for (size_t i = 0; i < b_values.size(); ++i) {
 norm_b += get_number(multiply_values(b_values[i], b_values[i]));
 }
 norm_b = std::sqrt(norm_b);
 // Вычисляем косинус угла
 double cos_angle = dot_product_ab / (norm_a * norm_b);
 // Ограничиваем значение косинуса в пределах [-1, 1], чтобы избежать 
ошибок из-за погрешностей вычислений
 cos_angle = std::min(1.0, std::max(-1.0, cos_angle));
 return Value(cos_angle);
}
Value _vector_mixed_product(Module_interface *, const VariableSetWrapper 
&args) {
 const Value &a_value = args[0].origin().value();
 const Value &b_value = args[1].origin().value();
 const Value &c_value = args[2].origin().value();
 assert(a_value.isArray() && b_value.isArray() && c_value.isArray() &&
 is_valid_array(*a_value.getArray()) && 
is_valid_array(*b_value.getArray()) && is_valid_array(*c_value.getArray()));
 const auto &a_array = a_value.getArray();
 const auto &b_array = b_value.getArray();
 const auto &c_array = c_value.getArray();
 const auto &a_values = a_array->values();
 const auto &b_values = b_array->values();
 const auto &c_values = c_array->values();
 // Проверяем, что все векторы трехмерные
 if (a_values.size() != 3 || b_values.size() != 3 || c_values.size() != 3) {
 throw bormental::DrBormental(
 "Module linalg vector_mixed_product",
 "All vectors must be three-dimensional"
 );
 }
 // Векторное произведение b и c
 double cross_product_b_c[3];
 cross_product_b_c[0] = get_number(b_values[1]) * get_number(c_values[2]) -
get_number(b_values[2]) * get_number(c_values[1]);
 cross_product_b_c[1] = get_number(b_values[2]) * get_number(c_values[0]) -
get_number(b_values[0]) * get_number(c_values[2]);
 cross_product_b_c[2] = get_number(b_values[0]) * get_number(c_values[1]) -
get_number(b_values[1]) * get_number(c_values[0]);
 // Скалярное произведение a и (b × c)
 double mixed_product = get_number(a_values[0]) * cross_product_b_c[0] +
 get_number(a_values[1]) * cross_product_b_c[1] +
 get_number(a_values[2]) * cross_product_b_c[2];
 return Value(mixed_product);
}
Value _matr_pow(Module_interface *, const VariableSetWrapper &args) {
 const Value &matrix_value = args[0].origin().value();
 const Value &exponent_value = args[1].origin().value();
 
 assert(matrix_value.isArray() && is_valid_array(*matrix_value.getArray()));
 assert(exponent_value.isInt());
 const auto &matrix_array = matrix_value.getArray();
 const auto &dimensions = matrix_array->dimensions();
 index_t exponent = exponent_value.getInt();
 if (dimensions.size() != 2 || dimensions[0] != dimensions[1]) {
 throw bormental::DrBormental(
 "Module linalg matr_pow",
 "The argument is not a square 2D array.");
 }
 index_t n = dimensions[0];
 const auto &values = matrix_array->values();
 std::vector<Value> result_matrix(n * n, Value(0));
 for (index_t i = 0; i < n; ++i) {
 result_matrix[i * n + i] = Value(1);
 }
 if (exponent == 0) {
 return std::make_shared<Array>(
 ValueType::Null,
 dimensions,
 result_matrix
 );
 }
 
 //if (exponent < 0) {
 // throw bormental::DrBormental(
 // "Module linalg matr_pow",
 // "Negative powers are not supported"
 // }
 
 std::vector<Value> base_matrix = values;
 while (exponent > 0) {
 if (exponent % 2 == 1) {
 std::vector<Value> temp_result(n*n, Value(0));
 for (index_t i=0; i<n; ++i) {
 for (index_t j=0; j < n; ++j) {
 Value sum(0);
 for (index_t k=0; k<n; ++k){
 sum = plus_values(sum, 
multiply_values(result_matrix[i*n+k], base_matrix[k * n+j]));
 }
 temp_result[i*n + j] = sum;
 }
 }
 result_matrix = temp_result;
 }
 
 std::vector<Value> temp_base(n*n, Value(0));
 for (index_t i = 0; i < n; ++i) {
 for (index_t j=0; j<n; ++j) {
 Value sum(0);
 for (index_t k = 0; k<n; ++k) {
 sum = plus_values(sum, 
multiply_values(base_matrix[i*n+k], base_matrix[k*n+j]));
 }
 temp_base[i*n+j] = sum;
 }
 }
 base_matrix = temp_base;
 
 exponent /= 2;
 }
 auto result = std::make_shared<Array>(
 ValueType::Null,
 dimensions,
 result_matrix
 );
 return result;
}
Value _matr_rank(Module_interface *, const VariableSetWrapper &args){
const Value &matrix_value = args[0].origin().value();
assert(matrix_value.isArray());
const auto &matrix_array = matrix_value.getArray();
const auto &dimensions = matrix_array->dimensions();
if (dimensions.size() != 2) {
throw bormental::DrBormental(
"Module linalg matr_rank",
"The argument must be a 2D matrix.");
}
index_t rows = dimensions[0];
index_t cols = dimensions[1];
auto matrix = matrix_array->values();
index_t rank = 0;
for (index_t col = 0; col < cols; ++col) {
index_t pivot_row = rank;
while (pivot_row < rows && 
real_is_zero<double>(get_number(matrix[pivot_row*cols+col]))) {
++pivot_row;
}
if (pivot_row == rows) {
continue;
}
if (pivot_row != rank) {
for (index_t j = 0; j < cols; ++j) {
std::swap(matrix[rank * cols + j], matrix[pivot_row * cols 
+ j]);
}
}
Value pivot = matrix[rank * cols + col];
if (!real_is_zero<double>(get_number(pivot))) {
for (index_t j = 0; j < cols; ++j) {
matrix[rank*cols+j] = 
get_number(matrix[rank*cols+j])/get_number(pivot);
}
}
for(index_t i = rank + 1; i < rows; ++i) {
Value factor = matrix[i*cols + col];
for (index_t j = 0; j < cols; ++j) {
Value mul = multiply_values(matrix[rank*cols+j],factor);
matrix[i*cols+j] = minus_values(matrix[i*cols+j], mul); 
// matrix[i*cols+j]
}
}
++rank;
}
return Value(int64_t(rank));
}
Value _solve_linear_system(Module_interface *, const VariableSetWrapper &args) 
{
 // Получаем матрицу коэффициентов и вектор правой части
 const Value &matrix_value = args[0].origin().value();
 const Value &rhs_value = args[1].origin().value();
 // Проверяем входные данные
 assert(matrix_value.isArray());
 assert(rhs_value.isArray());
 const auto &matrix_array = matrix_value.getArray();
 const auto &rhs_array = rhs_value.getArray();
 const auto &matrix_dimensions = matrix_array->dimensions();
 const auto &rhs_dimensions = rhs_array->dimensions();
 if (matrix_dimensions.size() != 2 || matrix_dimensions[0] != rhs_dimensions[0]) 
{
 throw bormental::DrBormental(
 "Module linalg solve_linear_system",
 "The coefficient matrix must be 2D and compatible with the RHS vector.");
 }
 index_t n = matrix_dimensions[0];
 index_t m = matrix_dimensions[1];
 if (n != m) {
 throw bormental::DrBormental(
 "Module linalg solve_linear_system",
 "The coefficient matrix must be square.");
 }
 const auto &matrix_values = matrix_array->values();
 const auto &rhs_values = rhs_array->values();
 // Приводим значения к стандартным типам
 std::vector<std::vector<double>> a(n, std::vector<double>(n));
 std::vector<double> b(n);
 for (index_t i = 0; i < n; ++i) {
 for (index_t j = 0; j < n; ++j) {
 a[i][j] = get_number(matrix_values[i * n + j]);
 }
 b[i] = get_number(rhs_values[i]);
 }
 // Прямой ход метода Гаусса
 for (index_t i = 0; i < n; ++i) {
 // Поиск ведущего элемента
 index_t max_row = i;
 for (index_t k = i + 1; k < n; ++k) {
 if (std::abs(a[k][i]) > std::abs(a[max_row][i])) {
 max_row = k;
 }
 }
 // Перестановка строк
 std::swap(a[i], a[max_row]);
 std::swap(b[i], b[max_row]);
 // Проверка на вырожденность
 if (std::abs(a[i][i]) < 1e-9) {
 throw bormental::DrBormental(
 "Module linalg solve_linear_system",
 "The system has no unique solution (matrix is singular).");
 }
 // Нормализация строки
 double diag_element = a[i][i];
 for (index_t j = 0; j < n; ++j) {
 a[i][j] /= diag_element;
 }
 b[i] /= diag_element;
 // Обнуление элементов столбца ниже диагонали
 for (index_t k = i + 1; k < n; ++k) {
 double factor = a[k][i];
 for (index_t j = 0; j < n; ++j) {
 a[k][j] -= factor * a[i][j];
 }
 b[k] -= factor * b[i];
 }
 }
 // Обратный ход метода Гаусса
 std::vector<double> x(n);
 for (int i = n - 1; i >= 0; --i) {
 x[i] = b[i];
 for (index_t j = i + 1; j < n; ++j) {
 x[i] -= a[i][j] * x[j];
 }
 }
 // Преобразуем результат в Value
 std::vector<Value> result_values;
 result_values.reserve(n);
 for (const auto &val : x) {
 result_values.emplace_back(Value(val));
 }
 return std::make_shared<Array>(
 ValueType::Null,
 std::vector<index_t>{index_t(n)},
 result_values);
}
Value _gram_matrix(Module_interface *, const VariableSetWrapper &args) {
 // Получение входной матрицы
 const Value &matrix_value = args[0].origin().value();
 // Проверка, что входные данные являются массивом (матрицей)
 assert(matrix_value.isArray());
 const auto &matrix_array = matrix_value.getArray();
 const auto &dimensions = matrix_array->dimensions();
 // Проверка, что матрица двумерная
 if (dimensions.size() != 2) {
 throw bormental::DrBormental(
 "Module linalg gram_matrix",
 "The input must be a 2D array (matrix).");
 }
 index_t n_rows = dimensions[0];
 index_t n_cols = dimensions[1];
 const auto &values = matrix_array->values();
 // Создание пустой матрицы Грама размером m x m
 std::vector<Value> gram_values(n_cols * n_cols, Value(0));
 // Вычисление элементов матрицы Грама
 for (index_t i = 0; i < n_cols; ++i) {
 for (index_t j = 0; j <= i; ++j) { // Используем симметрию
 Value scalar_product = Value(0);
 for (index_t k = 0; k < n_rows; ++k) {
 Value subtotal = multiply_values(values[k * n_cols + i],values[k * 
n_cols + j]);
 scalar_product = plus_values(scalar_product,subtotal);
 }
 gram_values[i * n_cols + j] = scalar_product;
 if (i != j) {
 gram_values[j * n_cols + i] = scalar_product; // Симметричный 
элемент
 }
 }
 }
 // Возвращаем матрицу Грама
 return std::make_shared<Array>(
 ValueType::Null,
 std::vector<index_t>{index_t(n_cols), index_t(n_cols)},
 std::move(gram_values)
 );
}
Value _vector_norm(Module_interface *, const VariableSetWrapper &args) {
 // Получение входного вектора
 const Value &vector_value = args[0].origin().value();
 const Value &norm_type_value = args[1].origin().value();
 // Проверка, что входные данные являются массивом (вектором)
 assert(vector_value.isArray());
 const auto &vector_array = vector_value.getArray();
 const auto &dimensions = vector_array->dimensions();
 // Проверка, что массив одномерный
 if (dimensions.size() != 1) {
 throw bormental::DrBormental(
 "Module linalg vector_norm",
 "The input must be a 1D array (vector).");
 }
 // Проверка типа нормы
 assert(norm_type_value.isInt());
 int norm_type = norm_type_value.getInt();
 if (norm_type < 0 || norm_type > 2) {
 throw bormental::DrBormental(
 "Module linalg vector_norm",
 "Invalid norm type. Must be 0 (L1), 1 (L2), or 2 (L∞).");
 }
 const auto &values = vector_array->values();
 index_t n_elements = dimensions[0];
 // Переменные для вычисления
 double result = 0.0;
 switch (norm_type) {
 case 0: { // Норма L^1
 for (index_t i = 0; i < n_elements; ++i) {
 result += std::abs(get_number(values[i]));
 }
 break;
 }
 case 1: { // Норма L^2
 double l2_norm_squared = 0.0;
 for (index_t i = 0; i < n_elements; ++i) {
 double abs_value = std::abs(get_number(values[i]));
 l2_norm_squared += abs_value * abs_value;
 }
 result = std::sqrt(l2_norm_squared);
 break;
 }
 case 2: { // Норма L^∞
 double linf_norm = 0.0;
 for (index_t i = 0; i < n_elements; ++i) {
 linf_norm = std::max(linf_norm, std::abs(get_number(values[i])));
 }
 result = linf_norm;
 break;
 }
 default: {
 throw bormental::DrBormental(
 "Module linalg vector_norm",
 "Unexpected norm type.");
 }
 }
 // Возвращаем результат как скалярное значение
 return Value(result);
}
/*
void qr_decomposition(const std::vector<Value> &matrix, std::vector<Value> 
&q_matrix, std::vector<Value> &r_matrix, index_t n) {
 q_matrix.assign(n * n, Value(0));
 r_matrix.assign(n * n, Value(0));
 // Реализация ортогонализации Грама-Шмидта для Q и R
 for (index_t k = 0; k < n; ++k) {
 // Копируем k-тый столбец
 for (index_t i = 0; i < n; ++i) {
 q_matrix[i * n + k] = matrix[i * n + k];
 }
 // Ортогонализация
 for (index_t j = 0; j < k; ++j) {
 Value dot_product(0);
 for (index_t i = 0; i < n; ++i) {
 Value subtotal = multiply_values(q_matrix[i * n + j], q_matrix[i * n + 
k]);
 dot_product = plus_values(dot_product, subtotal);
 }
 for (index_t i = 0; i < n; ++i) {
 Value subtotal = multiply_values(dot_product, q_matrix[i * n + j]);
 q_matrix[i * n + k] = minus_values(q_matrix[i * n + k],subtotal);
 }
 }
 // Нормализация
 Value norm(0);
 for (index_t i = 0; i < n; ++i) {
 Value subtotal = multiply_values(q_matrix[i * n + k] ,q_matrix[i * n + k]);
 norm = plus_values(norm, subtotal);
 }
 norm = sqrt(get_number(norm));
 for (index_t i = 0; i < n; ++i) {
 double subtotal = get_number(q_matrix[i * n + k]) / get_number(norm);
 q_matrix[i * n + k] = subtotal;
 }
 }
 // Вычисляем R
 for (index_t i = 0; i < n; ++i) {
 for (index_t j = i; j < n; ++j) {
 Value dot_product(0);
 for (index_t k = 0; k < n; ++k) {
 Value subtotal = multiply_values(q_matrix[k * n + i],matrix[k * n + 
j]);
 dot_product = plus_values(dot_product, subtotal);
 }
 r_matrix[i * n + j] = dot_product;
 }
 }
}
bool has_converged(const std::vector<Value> &matrix_a, const std::vector<Value> 
&matrix_b, double tolerance, index_t n) {
 for (index_t i = 0; i < n; ++i) {
 for (index_t j = 0; j < n; ++j) {
 if (abs(get_number(minus_values(matrix_a[i * n + j],matrix_b[i * n + j]))) 
> tolerance) {
 return false;
 }
 }
 }
 return true;
}
Value _eigenvalues(Module_interface *, const VariableSetWrapper &args) {
 const Value &matrix_value = args[0].origin().value();
 assert(matrix_value.isArray());
 const auto &matrix_array = matrix_value.getArray();
 const auto &dimensions = matrix_array->dimensions();
 // Проверка, что матрица квадратная
 if (dimensions.size() != 2 || dimensions[0] != dimensions[1]) {
 throw bormental::DrBormental(
 "Module linalg eigenvalues",
 "The input matrix must be square (n x n).");
 }
 index_t n = dimensions[0];
 std::vector<Value> current_matrix = matrix_array->values();
 // QR-итерации
 const double tolerance = 1e-9;
 const size_t max_iterations = 1000;
 size_t iteration = 0;
 while (iteration < max_iterations) {
 std::vector<Value> q_matrix, r_matrix;
 // Выполнить QR-разложение
 qr_decomposition(current_matrix, q_matrix, r_matrix, n);
 // Обновляем матрицу
 std::vector<Value> next_matrix(n * n, Value(0));
 for (index_t i = 0; i < n; ++i) {
 for (index_t j = 0; j < n; ++j) {
 Value sum(0);
 for (index_t k = 0; k < n; ++k) {
 Value subtotal = multiply_values(r_matrix[i * n + k],q_matrix[k * n + 
j]);
 sum = plus_values(sum, subtotal);
 }
 next_matrix[i * n + j] = sum;
 }
 }
 // Проверка сходимости
 if (has_converged(current_matrix, next_matrix, tolerance, n)) {
 break;
 }
 current_matrix = next_matrix;
 ++iteration;
 }
 if (iteration == max_iterations) {
 throw bormental::DrBormental(
 "Module linalg eigenvalues",
 "QR method did not converge within the maximum number of iterations.");
 }
 // Извлекаем собственные значения из диагонали
 std::vector<Value> eigenvalues(n, Value(0));
 for (index_t i = 0; i < n; ++i) {
 eigenvalues[i] = current_matrix[i * n + i];
 }
 //return std::make_shared<Array>(
 // ValueType::Null,
 // std::vector<index_t>{index_t(n)},
 // std::move(eigenvalues)
 //);
 return std::make_shared<Array>(
 ValueType::Null,
 std::vector<index_t>{index_t(n)},
 std::move(current_matrix)
 );
 
}
*/
Value _eigenvalues(Module_interface *, const VariableSetWrapper &args) {
 const Value &matrix_value = args[0].origin().value();
 
 if (!matrix_value.isArray()) {
 throw std::runtime_error("Input must be an array");
 }
 
 const auto &matrix_array = matrix_value.getArray();
 const auto &dimensions = matrix_array->dimensions();
 
 if (dimensions.size() != 2 || dimensions[0] != dimensions[1]) {
 throw std::runtime_error("Matrix must be square");
 }
 
 index_t n = dimensions[0];
 std::vector<Value> matrix = matrix_array->values();
 
 const double tolerance = 1e-9;
 const size_t max_iterations = 100;
 
 for (size_t iter = 0; iter < max_iterations; ++iter) {
 std::vector<Value> q_matrix(n * n, Value(0));
 std::vector<Value> r_matrix(n * n, Value(0));
 
 // Ортогонализация Грамма-Шмидта
 for (index_t k = 0; k < n; ++k) {
 // Копирование k-го столбца
 for (index_t i = 0; i < n; ++i) {
 q_matrix[i * n + k] = matrix[i * n + k];
 }
 
 // Вычитание проекций
 for (index_t j = 0; j < k; ++j) {
 Value proj_scalar = Value(0);
 for (index_t i = 0; i < n; ++i) {
 proj_scalar = plus_values(proj_scalar, 
 multiply_values(q_matrix[i * n + j], q_matrix[i * n + k]));
 }
 
 for (index_t i = 0; i < n; ++i) {
 q_matrix[i * n + k] = minus_values(
 q_matrix[i * n + k], 
 multiply_values(proj_scalar, q_matrix[i * n + j])
 );
 }
 }
 
 // Безопасная нормализация
 Value norm_squared = Value(0);
 for (index_t i = 0; i < n; ++i) {
 norm_squared = plus_values(norm_squared, 
 multiply_values(q_matrix[i * n + k], q_matrix[i * n + k]));
 }
 
 double norm = sqrt(std::max(get_number(norm_squared), 1e-10));
 
 for (index_t i = 0; i < n; ++i) {
 q_matrix[i * n + k] = Value(get_number(q_matrix[i * n + k]) / norm);
 }
 }
 
 // Вычисление R-матрицы
 for (index_t i = 0; i < n; ++i) {
 for (index_t j = i; j < n; ++j) {
 Value dot_product = Value(0);
 for (index_t k = 0; k < n; ++k) {
 dot_product = plus_values(dot_product, 
 multiply_values(q_matrix[k * n + i], matrix[k * n + j]));
 }
 r_matrix[i * n + j] = dot_product;
 }
 }
 
 // Обновление матрицы
 std::vector<Value> next_matrix(n * n, Value(0));
 for (index_t i = 0; i < n; ++i) {
 for (index_t j = 0; j < n; ++j) {
 Value sum = Value(0);
 for (index_t k = 0; k < n; ++k) {
 sum = plus_values(sum, 
 multiply_values(r_matrix[i * n + k], q_matrix[k * n + j]));
 }
 next_matrix[i * n + j] = sum;
 }
 }
 
 // Проверка сходимости с защитой от переполнения
 bool converged = true;
 for (index_t i = 0; i < n; ++i) {
 for (index_t j = 0; j < n; ++j) {
 double diff = std::abs(get_number(minus_values(matrix[i * n + j], 
next_matrix[i * n + j])));
 if (diff > tolerance && diff < 1e10) {
 converged = false;
 break;
 }
 }
 if (!converged) break;
 }
 
 if (converged) {
 matrix = next_matrix;
 break;
 }
 
 matrix = next_matrix;
 }
 
 // Извлечение собственных значений
 std::vector<Value> eigenvalues(n, Value(0));
 for (index_t i = 0; i < n; ++i) {
 eigenvalues[i] = matrix[i * n + i];
 }
 
 return std::make_shared<Array>(
 ValueType::Null,
 std::vector<index_t>{index_t(n)},
 std::move(eigenvalues)
 );
}
Value _pinv(Module_interface *, const VariableSetWrapper &args) {
 const Value &matrix_value = args[0].origin().value();
 
 if (!matrix_value.isArray()) {
 throw std::runtime_error("Input must be an array");
 }
 
 const auto &matrix_array = matrix_value.getArray();
 const auto &dimensions = matrix_array->dimensions();
 
 if (dimensions.size() != 2) {
 throw std::runtime_error("Input must be a 2D matrix");
 }
 
 index_t rows = dimensions[0], cols = dimensions[1];
 std::vector<Value> input_matrix = matrix_array->values();
 
 // Epsilon для численной стабильности
 const double epsilon = 1e-10;
 
 // Преобразование входной матрицы
 std::vector<std::vector<double>> A(rows, std::vector<double>(cols));
 for (index_t i = 0; i < rows; ++i) {
 for (index_t j = 0; j < cols; ++j) {
 A[i][j] = get_number(input_matrix[i * cols + j]);
 }
 }
 
 // Внутренняя функция транспонирования
 auto transpose = [](const std::vector<std::vector<double>>& matrix) {
 size_t rows = matrix.size(), cols = matrix[0].size();
 std::vector<std::vector<double>> result(cols, std::vector<double>(rows, 0.0));
 for (size_t i = 0; i < rows; ++i) {
 for (size_t j = 0; j < cols; ++j) {
 result[j][i] = matrix[i][j];
 }
 }
 return result;
 };
 
 // Внутренняя функция умножения матриц
 auto multiply = [](const std::vector<std::vector<double>>& A, 
 const std::vector<std::vector<double>>& B) {
 size_t rowsA = A.size(), colsA = A[0].size();
 size_t rowsB = B.size(), colsB = B[0].size();
 if (colsA != rowsB) {
 throw std::invalid_argument("Matrix dimensions do not match for 
multiplication.");
 }
 std::vector<std::vector<double>> result(rowsA, std::vector<double>(colsB, 
0.0));
 for (size_t i = 0; i < rowsA; ++i) {
 for (size_t j = 0; j < colsB; ++j) {
 for (size_t k = 0; k < colsA; ++k) {
 result[i][j] += A[i][k] * B[k][j];
 }
 }
 }
 return result;
 };
 
 // SVD методом Якоби с оптимизацией
 auto svd = [&epsilon](const std::vector<std::vector<double>>& A) 
 -> std::tuple<std::vector<std::vector<double>>, 
 std::vector<std::vector<double>>, 
 std::vector<std::vector<double>>> {
 size_t m = A.size(), n = A[0].size();
 std::vector<std::vector<double>> U = A;
 std::vector<std::vector<double>> V(n, std::vector<double>(n, 0.0));
 for (size_t i = 0; i < n; ++i) V[i][i] = 1.0;
 std::vector<std::vector<double>> S(m, std::vector<double>(n, 0.0));
 bool converge = false;
 int max_iterations = 1000; // Ограничение итераций для предотвращения 
бесконечного цикла
 
 while (!converge && max_iterations-- > 0) {
 converge = true;
 for (size_t p = 0; p < n; ++p) {
 for (size_t q = p + 1; q < n; ++q) {
 double alpha = 0.0, beta = 0.0, gamma = 0.0;
 for (size_t i = 0; i < m; ++i) {
 alpha += U[i][p] * U[i][p];
 beta += U[i][q] * U[i][q];
 gamma += U[i][p] * U[i][q];
 }
 if (std::abs(gamma) > epsilon) {
 converge = false;
 double zeta = (beta - alpha) / (2.0 * gamma);
 double t = std::copysign(1.0 / (std::abs(zeta) + std::sqrt(1 + zeta * 
zeta)), zeta);
 double c = 1.0 / std::sqrt(1 + t * t);
 double s = c * t;
 for (size_t i = 0; i < m; ++i) {
 double up = U[i][p], uq = U[i][q];
 U[i][p] = c * up - s * uq;
 U[i][q] = s * up + c * uq;
 }
 for (size_t i = 0; i < n; ++i) {
 double vp = V[i][p], vq = V[i][q];
 V[i][p] = c * vp - s * vq;
 V[i][q] = s * vp + c * vq;
 }
 }
 }
 }
 }
 // Нормализация
 for (size_t j = 0; j < n; ++j) {
 double sigma = 0.0;
 for (size_t i = 0; i < m; ++i) {
 sigma += U[i][j] * U[i][j];
 }
 S[j][j] = std::sqrt(std::max(sigma, epsilon));
 for (size_t i = 0; i < m; ++i) {
 U[i][j] /= S[j][j];
 }
 }
 return {U, S, V};
 };
 
 // Вычисление псевдообратной матрицы
 auto compute_pseudoinverse = [&](const std::vector<std::vector<double>>& A) 
{
 auto [U, S, V] = svd(A);
 
 size_t m = S.size(), n = S[0].size();
 std::vector<std::vector<double>> S_inv(n, std::vector<double>(m, 0.0));
 
 for (size_t i = 0; i < std::min(m, n); ++i) {
 if (S[i][i] > epsilon) {
 S_inv[i][i] = 1.0 / S[i][i];
 }
 }
 auto V_T = transpose(V);
 auto U_T = transpose(U);
 return multiply(multiply(V, S_inv), U_T);
 };
 
 // Вычисление псевдообратной матрицы
 std::vector<std::vector<double>> A_pseudo = compute_pseudoinverse(A);
 
 // Конвертация результата обратно в Value
 std::vector<Value> result_matrix;
 for (const auto& row : A_pseudo) {
 for (double val : row) {
 result_matrix.push_back(Value(val));
 }
 }
 
 // Возврат как многомерного массива
 return std::make_shared<Array>(
 ValueType::Null,
 std::vector<index_t>{index_t(A_pseudo.size()), 
index_t(A_pseudo[0].size())},
 std::move(result_matrix)
 );
}
}