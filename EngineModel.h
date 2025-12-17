#pragma once
#include "Types.h"
#include <cmath>
#include <stdexcept>

class EngineModel {
public:
    // ==========================================
    // ПАРАМЕТРИ СИСТЕМИ з рівняння (1)
    // 
    // Рівняння (1):
    // T·d⁴x/dt⁴ + (1 + r·T·k₂)·d³x/dt³ + T·k₁·k₂·k₃·d²x/dt² = 
    //     = k₁·T·d³F/dt³ + (k₁ + r·T·k₂)·d²F/dt²
    // ==========================================

    // Постійна часу
    static constexpr double T = 0.1;
    
    // Коефіцієнт зворотного зв'язку
    static constexpr double r = 1.5;
    
    // Коефіцієнти передачі ланок системи
    static constexpr double k1 = 2.0;
    static constexpr double k2 = 1.0;
    static constexpr double k3 = 0.5;

    // ==========================================
    // ФУНКЦІЯ ЗБУРЕННЯ F(t)
    // Використовуємо гладку функцію, яка має всі похідні
    // F(t) = A·exp(-α·t)·sin(ω·t)  для t ≥ 0
    // ==========================================
    
    static constexpr double A = 1.0;      // Амплітуда
    static constexpr double alpha = 0.5;  // Коефіцієнт загасання
    static constexpr double omega = 3.0;  // Частота коливань

    static double F(double t) {
        if (t < 0) return 0.0;
        return A * std::exp(-alpha * t) * std::sin(omega * t);
    }

    // ==========================================
    // ПОХІДНІ ЗБУРЕННЯ F(t)
    // F(t) = A·e^(-α·t)·sin(ω·t)
    //
    // F'(t) = A·e^(-α·t)·[ω·cos(ω·t) - α·sin(ω·t)]
    //
    // F''(t) = A·e^(-α·t)·[(α²-ω²)·sin(ω·t) - 2·α·ω·cos(ω·t)]
    //
    // F'''(t) = A·e^(-α·t)·[(3·α·ω²-α³)·sin(ω·t) + (ω³-3·α²·ω)·cos(ω·t)]
    // ==========================================
    
    static double F_first_derivative(double t) {
        if (t < 0) return 0.0;
        double exp_term = A * std::exp(-alpha * t);
        return exp_term * (omega * std::cos(omega * t) - alpha * std::sin(omega * t));
    }

    static double F_second_derivative(double t) {
        if (t < 0) return 0.0;
        double exp_term = A * std::exp(-alpha * t);
        return exp_term * ((alpha * alpha - omega * omega) * std::sin(omega * t) 
                          - 2.0 * alpha * omega * std::cos(omega * t));
    }

    static double F_third_derivative(double t) {
        if (t < 0) return 0.0;
        double exp_term = A * std::exp(-alpha * t);
        double alpha2 = alpha * alpha;
        double alpha3 = alpha2 * alpha;
        double omega2 = omega * omega;
        double omega3 = omega2 * omega;
        
        return exp_term * ((3.0 * alpha * omega2 - alpha3) * std::sin(omega * t) 
                          + (omega3 - 3.0 * alpha2 * omega) * std::cos(omega * t));
    }

    // ==========================================
    // СИСТЕМА ДИФЕРЕНЦІАЛЬНИХ РІВНЯНЬ
    // 
    // Зводимо рівняння 4-го порядку до системи 1-го порядку:
    // 
    // y₀ = x        →  dy₀/dt = y₁
    // y₁ = x'       →  dy₁/dt = y₂
    // y₂ = x''      →  dy₂/dt = y₃
    // y₃ = x'''     →  dy₃/dt = x⁽⁴⁾
    //
    // З рівняння (1) виражаємо x⁽⁴⁾:
    // T·x⁽⁴⁾ = k₁·T·F⁽³⁾ + (k₁ + r·T·k₂)·F⁽²⁾ 
    //          - (1 + r·T·k₂)·x⁽³⁾ - T·k₁·k₂·k₃·x⁽²⁾
    //
    // x⁽⁴⁾ = [k₁·T·F⁽³⁾ + (k₁ + r·T·k₂)·F⁽²⁾ 
    //         - (1 + r·T·k₂)·x⁽³⁾ - T·k₁·k₂·k₃·x⁽²⁾] / T
    // ==========================================
    static State computeDerivatives(double t, const State& current_state) {
        State dydt;

        // Розпакування вектора стану
        // double x         = current_state[0];  // не потрібно для обчислення
        // double x_prime   = current_state[1];  // не потрібно для обчислення
        double x_double  = current_state[2];  // x''  (друга похідна)
        double x_triple  = current_state[3];  // x''' (третя похідна)

        // --- Кінематичні зв'язки (зниження порядку) ---
        dydt[0] = current_state[1];  // dx/dt = x'
        dydt[1] = current_state[2];  // d(x')/dt = x''
        dydt[2] = current_state[3];  // d(x'')/dt = x'''

        // --- Динаміка (з рівняння 1) ---
        // Отримуємо похідні збурення
        double F_dd = F_second_derivative(t);   // d²F/dt²
        double F_ddd = F_third_derivative(t);   // d³F/dt³

        // Обчислюємо коефіцієнти
        double coef_x_triple = (1.0 + r * T * k2);     // Коефіцієнт при x'''
        double coef_x_double = T * k1 * k2 * k3;        // Коефіцієнт при x''
        double coef_F_dd = (k1 + r * T * k2);           // Коефіцієнт при F''
        double coef_F_ddd = k1 * T;                     // Коефіцієнт при F'''

        // Перевірка на сингулярність
        if (T == 0.0) {
            throw std::runtime_error("T cannot be zero (singular system).");
        }

        // Обчислюємо x⁽⁴⁾ з рівняння (1)
        double x_fourth = (coef_F_ddd * F_ddd + coef_F_dd * F_dd 
                          - coef_x_triple * x_triple 
                          - coef_x_double * x_double) / T;

        dydt[3] = x_fourth;  // d(x''')/dt = x⁽⁴⁾

        return dydt;
    }

    // ==========================================
    // ДОПОМІЖНІ ФУНКЦІЇ для аналізу
    // ==========================================
    
    // Обчислення позначень з рівняння (2) для діагностики
    // (якщо знадобиться для розширеного завдання)
    static double compute_C1() {
        return T * k1 * k2 * k3;
    }

    static double compute_C2() {
        return (1.0 + r * T * k2);
    }

    static double compute_C3() {
        return T;
    }

    // Вивід параметрів системи
    static void printParameters() {
        std::cout << "=== Параметри САР ===" << std::endl;
        std::cout << "T  = " << T << " (постійна часу)" << std::endl;
        std::cout << "r  = " << r << " (коефіцієнт зворотного зв'язку)" << std::endl;
        std::cout << "k1 = " << k1 << " (коефіцієнт передачі 1)" << std::endl;
        std::cout << "k2 = " << k2 << " (коефіцієнт передачі 2)" << std::endl;
        std::cout << "k3 = " << k3 << " (коефіцієнт передачі 3)" << std::endl;
        std::cout << std::endl;
        std::cout << "Коефіцієнти в позначеннях (2):" << std::endl;
        std::cout << "C1 = " << compute_C1() << std::endl;
        std::cout << "C2 = " << compute_C2() << std::endl;
        std::cout << "C3 = " << compute_C3() << std::endl;
        std::cout << "=====================" << std::endl;
    }
};
