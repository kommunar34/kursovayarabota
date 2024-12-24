#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>

/**
 * @brief Класс, представляющий математический вектор с различными операциями.
 *
 * @tparam T Тип элементов вектора (например, int, double).
 */
template <typename T>
class Vector {
private:
    std::vector<T> values; ///< Элементы вектора.

public:
    /**
     * @brief Конструктор, создающий вектор из заданных данных.
     *
     * @param data std::vector, содержащий элементы вектора.
     * @throws std::invalid_argument если вектор пустой.
     */
    explicit Vector(std::vector<T> data) : values(std::move(data)) {
        if (values.empty()) {
            throw std::invalid_argument("Vector cannot be empty.");
        }    
}

    /**
     * @brief Возвращает количество элементов в векторе.
     * @return Размер вектора.
     */
    size_t size() const {
        return values.size();
    }

    /**
     * @brief Получает значение по заданному индексу.
     *
     * @param index Индекс элемента.
     * @return Значение элемента по указанному индексу.
     * @throws std::out_of_range если индекс выходит за пределы вектора.
     */
    T get(size_t index) const {
        if (index >= values.size()) {
            throw std::out_of_range("Index out of range.");
        }
        return values[index];
    }

    /**
     * @brief Устанавливает значение элемента по заданному индексу.
     *
     * @param index Индекс элемента.
     * @param value Новое значение для элемента.
     * @throws std::out_of_range если индекс выходит за пределы вектора.
     */
    void set(size_t index, T value) {
        if (index >= values.size()) {
            throw std::out_of_range("Index out of range.");
        }
        values[index] = value;
    }

    /**
     * @brief Вычисляет скалярное произведение с другим вектором.
     *
     * @param other Другой вектор.
     * @return Скалярное произведение двух векторов.
     * @throws std::invalid_argument если размеры векторов не совпадают.
     */
    T dot(const Vector<T>& other) const {
        if (size() != other.size()) {
            throw std::invalid_argument("Vectors must be of the same size for dot product.");
        }
        T result = T();
        for (size_t i = 0; i < size(); ++i) {
            result += values[i] * other.values[i];
        }
        return result;
    }

    /**         

     * @brief Вычисляет векторное произведение с другим вектором.
     *
     * @param other Другой вектор.
     * @return Новый вектор, представляющий векторное произведение.
     * @throws std::invalid_argument если размер любого из векторов не равен 3.
     */
    Vector<T> cross(const Vector<T>& other) const {
        if (size() != 3 || other.size() != 3) {
            throw std::invalid_argument("Cross product is defined only for 3D vectors.");
        }
        return Vector<T>({
            values[1] * other.values[2] - values[2] * other.values[1],
            values[2] * other.values[0] - values[0] * other.values[2],
            values[0] * other.values[1] - values[1] * other.values[0]
        });
    }

    /**
     * @brief Вычисляет смешанное произведение (тройное скалярное произведение) с двумя другими векторами.
     *
     * @param v2 Второй вектор.
     * @param v3 Третий вектор.
     * @return Скалярный результат смешанного произведения.   
     */
    T mixed_product(const Vector<T>& v2, const Vector<T>& v3) const {
        return this->dot(v2.cross(v3));

    }

    /**
     * @brief Вычисляет евклидову норму (длину) вектора.
     *
     * @return Норма вектора.
     */
    T norm(int norm_type = 1) const {
        T sum = T();
        switch (norm_type) {
            case 0: { // Норма L^1
                for (const auto& value : values) sum += std::abs(value);
            break;
            }
            case 1: { // Норма L^2
                for (const auto& value : values) sum += value * value;
                sum = std::sqrt(sum);
            break;
            }
            case 2: { // Норма L^∞
                for (const auto& value : values) sum = std::max(sum, std::abs(value));
            break;
            }
            default: {
                throw std::invalid_argument("Unexpected norm type.");
                }
            }
        return sum;
    }

    /**
     * @brief Вычисляет косинус угла между этим вектором и другим вектором.
     *
     * @param other Другой вектор.
     * @return Косинус угла между двумя векторами.
     */
    T cosine_with(const Vector<T>& other) const {
        if (norm() == 0 || other.norm() == 0) {
            throw std::invalid_argument("Norms of vectors must be are different from zero");
        }
        return dot(other) / (norm() * other.norm());
    }
    /**        
     * @brief Вычисляет скалярную проекцию этого вектора на другой вектор.
     *
     * @param other Другой вектор.
     * @return Значение скалярной проекции.
     */
    T projection_onto(const Vector<T>& other) const {
        if (norm() == 0 || other.norm() == 0) {
            throw std::invalid_argument("Norms of vectors must be are different from zero");
        }
        return dot(other) / other.norm();
    }
};
