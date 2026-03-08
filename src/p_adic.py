"""
@file p_adic.py
@brief p-adische Zahlen – Theorie und Implementierung.
@description
    Implementiert p-adische Zahlen, die Metrik und grundlegende Operationen.
    p-adische Zahlen sind ein alternatives Zahlensystem, das auf einem
    anderen Abstandsbegriff basiert als die reellen Zahlen.

    Kernidee: Der p-adische Betrag |n|_p = p^(-v_p(n)) wobei v_p(n) die
    p-adische Bewertung (die höchste Potenz von p, die n teilt) ist.

    Statt "große Zahlen sind weit entfernt" gilt: "durch hohe p-Potenzen
    teilbare Zahlen sind nah bei 0". Beispiel mit p=2:
        |8|_2 = 2^{-3} = 1/8  (klein!)
        |3|_2 = 2^0 = 1       (normal)
        |1/4|_2 = 2^2 = 4     (groß!)

    Satz von Ostrowski: Alle nicht-trivialen absoluten Werte auf Q sind
    entweder der gewöhnliche Betrag |·|_∞ oder ein p-adischer Betrag |·|_p.

    Anwendungen:
        - Hensels Lemma (p-adisches Newton-Verfahren)
        - Lokal-globale Prinzipien (Hasse-Minkowski)
        - p-adische L-Funktionen
        - Kryptographie (p-adische Gitter)

@author Kurt Ingwer
@version 1.0
@since 2026-03-08
@lastModified 2026-03-08
"""

import math
from typing import Union


# ===========================================================================
# GRUNDLEGENDE P-ADISCHE FUNKTIONEN
# ===========================================================================

def p_adic_valuation(n: int, p: int) -> Union[int, float]:
    """
    Berechnet die p-adische Bewertung v_p(n).

    Definition:
        v_p(n) = max{k ∈ ℤ : p^k | n}
               = höchste Potenz von p, die n teilt

    Beispiele:
        v_2(12) = 2   (12 = 4 · 3 = 2² · 3)
        v_3(12) = 1   (12 = 4 · 3 = 2² · 3¹)
        v_5(12) = 0   (5 ∤ 12)
        v_2(0) = +∞   (0 ist durch jede p-Potenz teilbar)

    Eigenschaften:
        v_p(mn) = v_p(m) + v_p(n)     (Bewertung ist ein Homomorphismus)
        v_p(m+n) ≥ min(v_p(m), v_p(n)) (ultrametrische Ungleichung)

    @param n: Ganzzahl (kann 0 oder negativ sein)
    @param p: Primzahl p ≥ 2
    @return: p-adische Bewertung v_p(n), oder float('inf') für n = 0
    @raises ValueError: Wenn p < 2
    @lastModified: 2026-03-08
    """
    if p < 2:
        raise ValueError(f"p muss eine Primzahl ≥ 2 sein, erhalten: p = {p}")

    # Sonderfall: 0 ist durch jede p-Potenz teilbar → v_p(0) = +∞
    if n == 0:
        return float('inf')

    # Arbeite mit absolutem Wert (Bewertung ist symmetrisch: v_p(-n) = v_p(n))
    n = abs(n)

    # Zähle, wie oft p in n aufgeht
    valuation = 0
    while n % p == 0:
        valuation += 1
        n //= p

    return valuation


def p_adic_norm(n: int, p: int) -> float:
    """
    Berechnet den p-adischen Betrag |n|_p.

    Definition:
        |n|_p = p^{-v_p(n)}   für n ≠ 0
        |0|_p = 0

    Der p-adische Betrag ist eine nicht-archimedische Norm:
        |n|_p ≥ 0                          (Positivität)
        |n|_p = 0  ⟺  n = 0              (Definitheit)
        |mn|_p = |m|_p · |n|_p            (Multiplikativität)
        |m+n|_p ≤ max(|m|_p, |n|_p)      (ultrametrische Dreiecksungleichung)

    Die letzte Eigenschaft ist stärker als die übliche Dreiecksungleichung!

    @param n: Ganzzahl
    @param p: Primzahl p ≥ 2
    @return: p-adischer Betrag ∈ [0, 1] für ganze Zahlen
    @lastModified: 2026-03-08
    """
    if n == 0:
        return 0.0

    v = p_adic_valuation(n, p)

    if v == float('inf'):
        return 0.0  # Sicherheitshalber (sollte durch n=0-Check oben abgedeckt sein)

    return float(p) ** (-v)


def p_adic_distance(x: int, y: int, p: int) -> float:
    """
    Berechnet den p-adischen Abstand d_p(x, y).

    Definition:
        d_p(x, y) = |x - y|_p = p^{-v_p(x-y)}

    Der p-adische Abstand ist eine ultrametrische Distanz:
        d_p(x, z) ≤ max(d_p(x, y), d_p(y, z))  (ultrametrische Ungleichung)

    Intuition: Zwei Zahlen sind p-adisch nah, wenn ihre Differenz
    durch eine hohe Potenz von p teilbar ist.

    Beispiel mit p=5:
        d_5(0, 25) = |25|_5 = 5^{-2} = 0.04   (sehr nah!)
        d_5(0, 3) = |3|_5 = 1                  (weit weg)

    @param x: Erste ganze Zahl
    @param y: Zweite ganze Zahl
    @param p: Primzahl p ≥ 2
    @return: p-adischer Abstand d_p(x, y) ≥ 0
    @lastModified: 2026-03-08
    """
    return p_adic_norm(x - y, p)


# ===========================================================================
# P-ADISCHE ZAHLEN ALS KLASSE
# ===========================================================================

class PAdicNumber:
    """
    Repräsentiert eine p-adische Zahl als unendliche p-adische Expansion.

    Eine p-adische Zahl hat die Form:
        x = Σ_{k=v}^∞ a_k · p^k
    wobei:
        v = v_p(x) die p-adische Bewertung ist
        a_k ∈ {0, 1, ..., p-1} die p-adischen Ziffern sind
        a_v ≠ 0 (erste signifikante Ziffer)

    Wir speichern endlich viele Ziffern als Approximation.

    Beispiel (p=2, Darstellung von -1):
        -1 = 1 + 1·2 + 1·4 + 1·8 + ... = Σ_{k=0}^∞ 2^k
        In 2-adischer Darstellung: ...1111111 (unendlich viele Einsen)
        Als endliche Approximation: [1,1,1,...,1] (mod 2^precision)

    @author Kurt Ingwer
    @since 2026-03-08
    @lastModified 2026-03-08
    """

    def __init__(self, digits: list, p: int, valuation: int = 0):
        """
        Initialisiert eine p-adische Zahl.

        @param digits: Liste der p-adischen Ziffern [a_v, a_{v+1}, ..., a_{v+n-1}]
                      Jede Ziffer muss in {0, ..., p-1} liegen.
        @param p: Primzahl (Basis der p-adischen Zahlen)
        @param valuation: p-adische Bewertung v (erster Exponent), kann negativ sein
        @raises ValueError: Wenn p < 2 oder Ziffern ungültig
        @lastModified: 2026-03-08
        """
        if p < 2:
            raise ValueError(f"p muss ≥ 2 sein, erhalten: p = {p}")
        for i, d in enumerate(digits):
            if not (0 <= d < p):
                raise ValueError(f"Ziffer an Position {i} muss in [0, {p-1}] liegen, erhalten: {d}")

        self.digits = list(digits)   # Kopie der Ziffernliste
        self.p = p                   # Basis
        self.valuation = valuation   # Bewertung (erster Exponent)

    @classmethod
    def from_integer(cls, n: int, p: int, precision: int = 20) -> 'PAdicNumber':
        """
        Konvertiert eine ganze Zahl in ihre p-adische Darstellung.

        Algorithmus:
            Wiederhole: Ziffer = n mod p, n = n // p
            Dies liefert die p-adischen Ziffern von niedrigster zu höchster Ordnung.

        Für negative Zahlen: Verwende das p-Komplement
            -n ≡ p^k - n (mod p^k) für hinreichend großes k

        @param n: Zu konvertierende ganze Zahl
        @param p: Primzahl
        @param precision: Anzahl der zu berechnenden Ziffern
        @return: PAdicNumber-Darstellung von n
        @lastModified: 2026-03-08
        """
        if n == 0:
            # Null hat keine signifikante Ziffer (alle Ziffern 0)
            return cls([0] * precision, p, valuation=0)

        # p-adische Bewertung bestimmen (Position der ersten Nicht-Null-Ziffer)
        val = p_adic_valuation(n, p)
        if val == float('inf'):
            val = 0  # Null-Fall (sollte oben abgedeckt sein)

        # Für negative Zahlen: konvertiere via p^k-Komplement
        # -n mod p^precision ergibt die p-adische Darstellung von -n
        if n < 0:
            # Berechne n mod p^precision (Python mod ist immer positiv)
            n_mod = n % (p ** precision)
        else:
            n_mod = n

        # Ziffernextraktion: n = a_0 + a_1·p + a_2·p² + ...
        digits = []
        temp = n_mod
        for _ in range(precision):
            digits.append(temp % p)
            temp //= p

        # Bewertung ist Position der ersten Nicht-Null-Ziffer
        actual_val = 0
        for i, d in enumerate(digits):
            if d != 0:
                actual_val = i
                break

        return cls(digits, p, valuation=actual_val)

    @classmethod
    def from_fraction(cls, num: int, den: int, p: int, precision: int = 20) -> 'PAdicNumber':
        """
        Konvertiert einen Bruch num/den in p-adische Darstellung.

        Methode: Berechne das modulare Inverse von den modulo p^precision,
        dann multipliziere mit num.

        Formel:
            num/den (mod p^k) = num · den^{-1} (mod p^k)

        Dies funktioniert, falls gcd(den, p) = 1 (d.h. p ∤ den).
        Falls p | den, hat der Bruch eine negative Bewertung (p-adischer Pol).

        @param num: Zähler
        @param den: Nenner (≠ 0)
        @param p: Primzahl
        @param precision: Anzahl der Ziffern
        @return: PAdicNumber-Darstellung von num/den
        @raises ValueError: Wenn den = 0 oder gcd(den, p) > 1 für einfache Fälle
        @lastModified: 2026-03-08
        """
        if den == 0:
            raise ValueError("Nenner darf nicht 0 sein")

        # Bestimme p-adische Bewertung des Nenners
        val_den = p_adic_valuation(den, p)
        val_num = p_adic_valuation(num, p)

        # Gesamtbewertung des Bruchs: v_p(num/den) = v_p(num) - v_p(den)
        val = val_num - (val_den if val_den != float('inf') else 0)
        if val == float('inf'):
            val = 0

        # Normalisiere: Teile num und den durch p^val_den (falls val_den endlich)
        if val_den != float('inf') and val_den > 0:
            p_pow = p ** int(val_den)
            while den % p_pow == 0:
                den //= p
        if val_num != float('inf') and val_num > 0:
            p_pow_num = p ** int(val_num)
            while num % p_pow_num == 0:
                num //= p

        # Jetzt: gcd(den_normalized, p) = 1, berechne Inverses mod p^precision
        mod = p ** precision
        try:
            den_inv = pow(den % mod, -1, mod)
        except (ValueError, Exception):
            # Falls kein Inverses existiert (p | den), approximiere
            den_inv = 1

        # num/den mod p^precision
        value = (num * den_inv) % mod

        # Extrahiere Ziffern
        digits = []
        temp = value
        for _ in range(precision):
            digits.append(temp % p)
            temp //= p

        # Tatsächliche Bewertung (erste Nicht-Null-Ziffer)
        actual_val = int(val) if val != float('inf') else 0

        return cls(digits, p, valuation=actual_val)

    def __add__(self, other: 'PAdicNumber') -> 'PAdicNumber':
        """
        Addiert zwei p-adische Zahlen mit Übertragsbehandlung.

        Addition erfolgt ziffer-für-ziffer wie bei normaler Ganzzahladdition,
        aber die Übertragspropagation geht "nach oben" (zu höheren p-Potenzen).

        @param other: Zweite p-adische Zahl (muss gleiches p haben)
        @return: Summe als PAdicNumber
        @raises ValueError: Wenn die Primbasen verschieden sind
        @lastModified: 2026-03-08
        """
        if self.p != other.p:
            raise ValueError(f"Primbasen müssen gleich sein: {self.p} ≠ {other.p}")

        p = self.p
        # Länge der Ausgabe: Maximum der beiden Zifernanzahlen + 1 (für Übertrag)
        n = max(len(self.digits), len(other.digits)) + 1

        # Fülle kürzere Liste mit Nullen auf
        d1 = self.digits + [0] * (n - len(self.digits))
        d2 = other.digits + [0] * (n - len(other.digits))

        # Ziffer-für-Ziffer-Addition mit Übertrag
        result_digits = []
        carry = 0
        for i in range(n):
            total = d1[i] + d2[i] + carry
            result_digits.append(total % p)
            carry = total // p

        # Bewertung der Summe
        new_val = min(self.valuation, other.valuation)

        return PAdicNumber(result_digits, p, valuation=new_val)

    def __mul__(self, other: 'PAdicNumber') -> 'PAdicNumber':
        """
        Multipliziert zwei p-adische Zahlen.

        Multiplikation erfolgt wie polynomiale Multiplikation mit anschließender
        Reduktion modulo p (Übertragspropagation).

        @param other: Zweite p-adische Zahl
        @return: Produkt als PAdicNumber
        @raises ValueError: Wenn die Primbasen verschieden sind
        @lastModified: 2026-03-08
        """
        if self.p != other.p:
            raise ValueError(f"Primbasen müssen gleich sein: {self.p} ≠ {other.p}")

        p = self.p
        n1 = len(self.digits)
        n2 = len(other.digits)
        n = n1 + n2  # Maximale Länge des Produkts

        # Konvolution der Ziffernlisten (wie Polynommultiplikation)
        raw = [0] * n
        for i in range(n1):
            for j in range(n2):
                raw[i + j] += self.digits[i] * other.digits[j]

        # Übertragspropagation
        result_digits = []
        carry = 0
        for i in range(n):
            total = raw[i] + carry
            result_digits.append(total % p)
            carry = total // p

        # Bewertung des Produkts: v_p(xy) = v_p(x) + v_p(y)
        new_val = self.valuation + other.valuation

        return PAdicNumber(result_digits, p, valuation=new_val)

    def to_integer_approx(self, terms: int = 10) -> int:
        """
        Konvertiert die p-adische Zahl zurück zu einer ganzen Zahl (Approximation).

        Berechnet: Σ_{k=0}^{terms-1} digits[k] · p^k (mod p^terms)

        Dies ergibt eine ganze Zahl im Bereich [0, p^terms - 1], die eine
        Approximation der p-adischen Zahl modulo p^terms ist.

        @param terms: Anzahl der zu verwendenden Ziffern (Genauigkeit)
        @return: Approximation als ganze Zahl im Bereich [0, p^terms - 1]
        @lastModified: 2026-03-08
        """
        result = 0
        p_power = 1  # p^k, beginnt bei p^0 = 1

        for k in range(min(terms, len(self.digits))):
            result += self.digits[k] * p_power
            p_power *= self.p

        return result

    def norm(self) -> float:
        """
        Berechnet die p-adische Norm dieser Zahl.

        Die digits-Liste enthält die p-adischen Ziffern von der niedrigsten
        zur höchsten Ordnung: digits[i] ist der Koeffizient von p^i.
        Die p-adische Bewertung v_p(x) ist die Position der ersten Nicht-Null-Ziffer.
        Damit gilt: |x|_p = p^{-v_p(x)}.

        @return: p-adische Norm p^{-v_p(x)}
        @lastModified: 2026-03-08
        """
        # Finde erste Nicht-Null-Ziffer → gibt die p-adische Bewertung
        for i, d in enumerate(self.digits):
            if d != 0:
                return float(self.p) ** (-i)
        # Alle Ziffern sind 0: Norm = 0
        return 0.0

    def __repr__(self) -> str:
        digits_str = ', '.join(str(d) for d in self.digits[:10])
        if len(self.digits) > 10:
            digits_str += ', ...'
        return f"PAdicNumber([{digits_str}], p={self.p}, val={self.valuation})"

    def __str__(self) -> str:
        # Zeige die Ziffern von höchster zu niedrigster Ordnung (mathematische Konvention)
        digits_display = self.digits[:15]
        if len(self.digits) > 15:
            digits_str = '...' + ''.join(str(d) for d in reversed(digits_display))
        else:
            digits_str = ''.join(str(d) for d in reversed(digits_display))
        return f"({digits_str})_{self.p}"


# ===========================================================================
# HENSELS LEMMA
# ===========================================================================

def hensel_lift(f_coeffs: list, p: int, initial_root: int, n_lifts: int = 5) -> int:
    """
    Hebt eine Wurzel von f modulo p auf modulo p^n via Hensels Lemma.

    Hensels Lemma (p-adisches Newton-Verfahren):
        Gegeben: f(a) ≡ 0 (mod p) und f'(a) ≢ 0 (mod p)
        Dann: Es gibt eindeutig ein b ≡ a (mod p) mit f(b) ≡ 0 (mod p²)

    Allgemeine Liftformel:
        a_{k+1} = a_k - f(a_k) / f'(a_k)   (in p-adischer Arithmetik)
               = a_k - f(a_k) · [f'(a_k)]^{-1}  (mod p^{k+1})

    Anwendungsbeispiel:
        Finde x mit x² ≡ -1 (mod 5^n):
        Startwurzel: x₀ = 2 (da 2² = 4 ≡ -1 mod 5)
        Lift: x₁ ≡ 7 (mod 25), x₂ ≡ 57 (mod 125), ...

    @param f_coeffs: Koeffizienten des Polynoms [a_n, ..., a_1, a_0]
                    (höchster Grad zuerst, wie numpy.polyval-Konvention)
                    Beispiel: x² + 1 → [1, 0, 1]
    @param p: Primzahl
    @param initial_root: Startapproximation (Wurzel mod p)
    @param n_lifts: Anzahl der Lift-Schritte (Genauigkeit mod p^{n_lifts+1})
    @return: Gehobene Wurzel modulo p^{n_lifts+1}
    @raises ValueError: Wenn f'(initial_root) ≡ 0 (mod p) (Liftbedingung verletzt)
    @lastModified: 2026-03-08
    """
    def poly_eval(coeffs: list, x: int, mod: int) -> int:
        """Wertet das Polynom an x aus (Horner-Schema, mod m)."""
        result = 0
        for c in coeffs:
            result = (result * x + c) % mod
        return result

    def poly_deriv_eval(coeffs: list, x: int, mod: int) -> int:
        """Wertet die formale Ableitung aus (Horner-Schema, mod m)."""
        n = len(coeffs) - 1  # Grad des Polynoms
        result = 0
        for i, c in enumerate(coeffs[:-1]):  # Konstante hat Ableitung 0
            # Koeffizient der Ableitung: (n-i) · coeffs[i] für Term x^{n-i}
            deriv_coeff = (n - i) * c
            result = (result * x + deriv_coeff) % mod
        return result

    # Prüfe Startbedingung: f(a) ≡ 0 (mod p)
    f_at_a = poly_eval(f_coeffs, initial_root, p)
    if f_at_a % p != 0:
        raise ValueError(
            f"Startbedingung verletzt: f({initial_root}) = {f_at_a} ≢ 0 (mod {p})"
        )

    # Prüfe Liftbedingung: f'(a) ≢ 0 (mod p)
    df_at_a = poly_deriv_eval(f_coeffs, initial_root, p)
    if df_at_a % p == 0:
        raise ValueError(
            f"Liftbedingung verletzt: f'({initial_root}) = {df_at_a} ≡ 0 (mod {p}). "
            f"Hensel-Lift nicht möglich (singuläre Wurzel)."
        )

    # Iterativer Lift: a_{k+1} = a_k - f(a_k) · f'(a_k)^{-1} (mod p^{k+1})
    current = initial_root
    current_mod = p

    for step in range(n_lifts):
        # Nächstes Modulus
        next_mod = current_mod * p

        # f(current) und f'(current) modulo next_mod berechnen
        f_val = poly_eval(f_coeffs, current, next_mod)
        df_val = poly_deriv_eval(f_coeffs, current, next_mod)

        # Modulares Inverses von f'(current) mod p (nach Liftbedingung existiert es)
        df_mod_p = df_val % p
        df_inv = pow(df_mod_p, -1, p)

        # Newton-Schritt: a_new = a - f(a) · f'(a)^{-1}  (mod p^{k+1})
        correction = (f_val * df_inv) % next_mod
        current = (current - correction) % next_mod

        current_mod = next_mod

    return current


# ===========================================================================
# P-ADISCHE ANALYTISCHE FUNKTIONEN
# ===========================================================================

def p_adic_exp(x_digits: list, p: int, precision: int = 10) -> list:
    """
    Berechnet die p-adische Exponentialfunktion e_p(x).

    Definition (via Potenzreihe):
        exp_p(x) = Σ_{n=0}^∞ x^n / n!

    Konvergenz: Die Reihe konvergiert p-adisch genau dann, wenn
        |x|_p < p^{-1/(p-1)}
    d.h. v_p(x) ≥ 1 für p ungerade, v_2(x) ≥ 2 für p = 2.

    Wir berechnen die Partialsumme bis zum angegebenen Grad.
    Die Koeffizienten werden als ganze Zahlen modulo p^precision gerechnet.

    @param x_digits: p-adische Ziffern von x (als Liste [a_0, a_1, ...])
    @param p: Primzahl
    @param precision: Anzahl der p-adischen Ziffern in der Ausgabe
    @return: Liste der p-adischen Ziffern von exp_p(x) (mod p^precision)
    @lastModified: 2026-03-08
    """
    mod = p ** precision

    # Konvertiere Ziffernliste zu ganzer Zahl (mod p^precision)
    x_int = 0
    for i, d in enumerate(x_digits[:precision]):
        x_int += d * (p ** i)
    x_int = x_int % mod

    # Partialsumme: Σ x^n / n! (mod p^precision)
    # Wir müssen sicherstellen, dass n! ein Inverses hat (also gcd(n!, p^k) = 1 für kleine n)
    result = 0
    x_power = 1    # x^n
    factorial = 1  # n!

    for n in range(precision * 3):  # Genug Terme für Konvergenz
        # n! muss zu p koprim sein für die ganzzahlige Inverskalkulation
        # Da p prim: p ∤ n! genau wenn n < p
        try:
            factorial_inv = pow(factorial % mod, -1, mod)
        except (ValueError, Exception):
            # factorial nicht invertierbar mod p^precision (n ≥ p)
            break

        term = (x_power * factorial_inv) % mod
        result = (result + term) % mod

        # Nächste Iteration
        x_power = (x_power * x_int) % mod
        factorial = factorial * (n + 1)
        if factorial > mod ** 2:  # Begrenze Wachstum
            factorial = factorial % mod

    # Extrahiere Ziffern
    digits = []
    temp = result
    for _ in range(precision):
        digits.append(temp % p)
        temp //= p

    return digits


def p_adic_log(x_digits: list, p: int, precision: int = 10) -> list:
    """
    Berechnet den p-adischen Logarithmus log_p(x).

    Definition (via Potenzreihe):
        log_p(x) = log_p(1 + (x-1)) = Σ_{n=1}^∞ (-1)^{n+1} (x-1)^n / n

    Konvergenz: Konvergiert p-adisch genau dann, wenn |x - 1|_p < 1,
    d.h. x ≡ 1 (mod p).

    @param x_digits: p-adische Ziffern von x (als Liste [a_0, a_1, ...])
                    WICHTIG: x muss ≡ 1 (mod p) sein für Konvergenz
    @param p: Primzahl
    @param precision: Anzahl der p-adischen Ziffern in der Ausgabe
    @return: Liste der p-adischen Ziffern von log_p(x) (mod p^precision)
    @lastModified: 2026-03-08
    """
    mod = p ** precision

    # Konvertiere x zu ganzer Zahl
    x_int = 0
    for i, d in enumerate(x_digits[:precision]):
        x_int += d * (p ** i)
    x_int = x_int % mod

    # y = x - 1 (muss |y|_p < 1 haben, d.h. p | y)
    y = (x_int - 1) % mod

    # Partialsumme: log(1+y) = Σ_{n=1}^∞ (-1)^{n+1} y^n / n
    result = 0
    y_power = y  # y^n, beginnt bei y^1

    for n in range(1, precision * 3 + 1):
        # n muss zu p koprim sein für Inverskalkulation (oder n < p)
        n_mod = n % mod
        if n_mod == 0:
            # n ist durch p^precision teilbar: Term verschwindet in mod-Arithmetik
            y_power = (y_power * y) % mod
            continue

        try:
            n_inv = pow(n_mod, -1, mod)
        except (ValueError, Exception):
            # Kein Inverses: überspringe Term
            y_power = (y_power * y) % mod
            continue

        sign = 1 if n % 2 == 1 else -1
        term = (sign * y_power * n_inv) % mod
        result = (result + term) % mod

        y_power = (y_power * y) % mod

    # Extrahiere Ziffern
    digits = []
    temp = result % mod
    for _ in range(precision):
        digits.append(temp % p)
        temp //= p

    return digits


# ===========================================================================
# SATZ VON OSTROWSKI (DEMONSTRATION)
# ===========================================================================

def ostrowski_theorem_demo(n: int) -> dict:
    """
    Demonstriert den Satz von Ostrowski für eine ganze Zahl n.

    Satz von Ostrowski (1916):
        Jeder nicht-triviale absolute Wert auf ℚ ist äquivalent zu einem der:
            |·|_∞  (gewöhnlicher absoluter Wert / archimedische Norm)
            |·|_p  für eine Primzahl p (p-adische Norm, nicht-archimedisch)

    Wichtige Folgerung – die Produktformel:
        |n|_∞ · Π_{p prim} |n|_p = 1   für alle n ∈ ℤ, n ≠ 0

    Diese "Produktformel" ist fundamental in der Zahlentheorie und verbindet
    alle lokalen Informationen (p-adische Beträge) mit der globalen Information
    (gewöhnlicher Betrag).

    Beispiel n = 12 = 2² · 3:
        |12|_∞ = 12
        |12|_2 = 2^{-2} = 1/4
        |12|_3 = 3^{-1} = 1/3
        |12|_5 = 5^0 = 1
        |12|_7 = 7^0 = 1
        Produkt: 12 · (1/4) · (1/3) · 1 · 1 · ... = 1 ✓

    @param n: Zu analysierende ganze Zahl (≠ 0)
    @return: Dictionary mit allen absoluten Werten und dem Produktnachweis
    @raises ValueError: Wenn n = 0
    @lastModified: 2026-03-08
    """
    if n == 0:
        raise ValueError("n muss ≠ 0 sein (|0|_p = 0 für alle p, Produktformel gilt nicht)")

    # Gewöhnlicher absoluter Wert
    abs_inf = abs(n)

    # p-adische Beträge für kleine Primzahlen
    small_primes = [2, 3, 5, 7, 11, 13, 17, 19]
    p_adic_values = {}

    for p in small_primes:
        p_adic_values[p] = p_adic_norm(n, p)

    # Produktformel: |n|_∞ · Π_p |n|_p sollte = 1 sein
    # (Über alle Primzahlen, aber für ganze n sind fast alle |n|_p = 1)
    product = abs_inf
    for p, val in p_adic_values.items():
        if val > 0:
            product *= val

    # Alle Primteiler von n finden (für vollständige Produktformel)
    prime_factors = {}
    temp = abs(n)
    for p in range(2, int(math.sqrt(abs(n))) + 2):
        if temp % p == 0:
            prime_factors[p] = 0
            while temp % p == 0:
                prime_factors[p] += 1
                temp //= p
    if temp > 1:
        prime_factors[temp] = 1

    # Produktformel über ALLE Primteiler
    exact_product = float(abs_inf)
    for p, exp in prime_factors.items():
        # |n|_p = p^{-exp}
        exact_product *= float(p) ** (-exp)

    return {
        'n': n,
        'abs_infinity': abs_inf,
        'p_adic_norms': p_adic_values,
        'prime_factorization': prime_factors,
        'product_formula_approx': product,      # |n|_∞ · Π_p |n|_p (kleine p)
        'product_formula_exact': exact_product, # Exaktes Produkt über alle Primteiler
        'product_equals_one': abs(exact_product - 1.0) < 1e-10,
        'ostrowski_explanation': (
            f"|{n}|_∞ = {abs_inf}, "
            + ", ".join(f"|{n}|_{p} = {v:.6f}" for p, v in p_adic_values.items() if v != 1.0)
            + f" → Produkt = {exact_product:.10f}"
        )
    }
