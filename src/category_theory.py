"""
@file category_theory.py
@brief Kategorientheorie: Objekte, Morphismen, Funktoren, natürliche Transformationen.
@description
    Dieses Modul implementiert die grundlegenden Strukturen der Kategorientheorie:

    Klassen:
    - Object               – Objekt einer Kategorie
    - Morphism             – Morphismus zwischen zwei Objekten
    - Category             – Abstrakte Kategorie mit Objekten und Morphismen
    - Funktor              – Kovarianter Funktor zwischen zwei Kategorien
    - NaturalTransformation – Natürliche Transformation zwischen zwei Funktoren

    Freie Funktionen:
    - set_category()         – Erstellt die Set-Kategorie für gegebene Mengen
    - group_category()       – Erstellt Grp-Kategorie aus Group-Objekten
    - opposite_category()    – Erstellt die duale Kategorie (alle Pfeile umgekehrt)
    - product_category()     – Produktkategorie C₁ × C₂
    - is_isomorphism()       – Prüft ob f ein Isomorphismus ist
    - is_epimorphism()       – Prüft ob f ein Epimorphismus ist
    - is_monomorphism()      – Prüft ob f ein Monomorphismus ist
    - hom_set()              – Gibt alle Morphismen von A nach B zurück
    - yoneda_functor()       – Yoneda-Einbettung h_A = Hom(A, -)

    Mathematische Grundlagen:
    - Eine Kategorie C besteht aus: Ob(C), Mor(C), Komposition ∘, Identitäten id_A
    - Assoziativität: h ∘ (g ∘ f) = (h ∘ g) ∘ f
    - Einheitsgesetze: id_B ∘ f = f = f ∘ id_A für f: A → B
    - Funktor F: C → D erhält Identitäten und Kompositionen
    - Nat. Transformation α: F ⇒ G mit Naturquadraten: α_B ∘ F(f) = G(f) ∘ α_A

@author Kurt Ingwer
@lastModified 2026-03-10
"""

from __future__ import annotations
from typing import Any, Callable, Optional
from itertools import product as iterproduct


# =============================================================================
# KLASSE: Object
# =============================================================================

class Object:
    """
    @brief Objekt einer Kategorie.
    @description
        Repräsentiert ein Objekt in einer Kategorie. Objekte werden durch ihren
        Namen identifiziert. In der Kategorientheorie sind Objekte abstrakte
        Entitäten – ihre interne Struktur spielt für die Kategorie keine Rolle.

        Mathematisch: Ob(C) ist die Klasse aller Objekte der Kategorie C.
    @lastModified 2026-03-10
    """

    def __init__(self, name: str, data: Any = None):
        """
        @brief Erstellt ein Kategorie-Objekt.
        @param name  Eindeutiger Bezeichner des Objekts
        @param data  Optionale interne Daten (z.B. eine Menge oder Gruppe)
        @lastModified 2026-03-10
        """
        # Name des Objekts (Bezeichner)
        self.name = name
        # Interne Daten (optional, z.B. eine Menge für Set-Kategorie)
        self.data = data

    def __repr__(self) -> str:
        """Lesbare Darstellung des Objekts."""
        return f"Object({self.name!r})"

    def __eq__(self, other: object) -> bool:
        """Gleichheit basiert auf dem Namen."""
        if isinstance(other, Object):
            return self.name == other.name
        return NotImplemented

    def __hash__(self) -> int:
        """Hash basiert auf dem Namen."""
        return hash(self.name)


# =============================================================================
# KLASSE: Morphism
# =============================================================================

class Morphism:
    """
    @brief Morphismus zwischen zwei Objekten einer Kategorie.
    @description
        Ein Morphismus f: A → B ist ein gerichteter Pfeil von einem
        Quellobjekt (source) zu einem Zielobjekt (target).

        Mathematisch:
        - f ∈ Hom(A, B) = Mor_C(A, B)
        - Komposition: (g ∘ f): A → C für f: A → B, g: B → C
        - Identität: id_A ∈ Hom(A, A) mit id_A ∘ f = f und g ∘ id_A = g
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        source: Object,
        target: Object,
        name: str,
        func: Optional[Callable] = None
    ):
        """
        @brief Erstellt einen Morphismus.
        @param source Quellobjekt A
        @param target Zielobjekt B
        @param name   Bezeichner des Morphismus (z.B. "f", "id_A")
        @param func   Optionale Python-Funktion, die den Morphismus repräsentiert
        @lastModified 2026-03-10
        """
        # Quell- und Zielobjekt
        self.source = source
        self.target = target
        # Name des Morphismus
        self.name = name
        # Zugehörige Python-Funktion (optional für konkrete Kategorien)
        self.func = func

    def compose(self, other: 'Morphism') -> 'Morphism':
        """
        @brief Komposition: self ∘ other (other zuerst, dann self).
        @description
            Berechnet die Komposition g ∘ f, wobei other = f: A → B
            und self = g: B → C. Das Ergebnis ist h = g ∘ f: A → C.

            Mathematisch: (g ∘ f)(x) = g(f(x))

            Vorbedingung: other.target == self.source
        @param other Morphismus f: A → B
        @return Neuer Morphismus g ∘ f: A → C
        @raises ValueError falls Quell/Ziel nicht übereinstimmen
        @lastModified 2026-03-10
        """
        # Prüfe Kompatibilität: Ziel von other muss Quelle von self sein
        if other.target != self.source:
            raise ValueError(
                f"Komposition nicht definiert: {other.name}.target={other.target.name} "
                f"≠ {self.name}.source={self.source.name}"
            )

        # Erstelle kombinierte Funktion, falls beide Funktionen vorhanden
        composed_func = None
        if self.func is not None and other.func is not None:
            # Schließe über lokale Variablen ab, um Closures korrekt zu erstellen
            f_func = other.func
            g_func = self.func
            composed_func = lambda x: g_func(f_func(x))

        return Morphism(
            source=other.source,
            target=self.target,
            name=f"({self.name}∘{other.name})",
            func=composed_func
        )

    def __repr__(self) -> str:
        """Lesbare Darstellung des Morphismus."""
        return f"Morphism({self.name!r}: {self.source.name} → {self.target.name})"

    def __eq__(self, other: object) -> bool:
        """Gleichheit basiert auf Name, Quelle und Ziel."""
        if isinstance(other, Morphism):
            return (self.name == other.name and
                    self.source == other.source and
                    self.target == other.target)
        return NotImplemented

    def __hash__(self) -> int:
        """Hash basiert auf Name, Quelle und Ziel."""
        return hash((self.name, self.source.name, self.target.name))


# =============================================================================
# KLASSE: Category
# =============================================================================

class Category:
    """
    @brief Abstrakte Kategorie mit Objekten und Morphismen.
    @description
        Eine Kategorie C besteht aus:
        1. Einer Klasse von Objekten Ob(C)
        2. Für je zwei Objekte A, B einer Menge Hom(A, B) von Morphismen
        3. Kompositionsoperationen ∘: Hom(B,C) × Hom(A,B) → Hom(A,C)
        4. Für jedes Objekt A einem Identitätsmorphismus id_A ∈ Hom(A,A)

        Axiome:
        - Assoziativität: (h ∘ g) ∘ f = h ∘ (g ∘ f)
        - Linke Einheit: id_B ∘ f = f
        - Rechte Einheit: f ∘ id_A = f
    @lastModified 2026-03-10
    """

    def __init__(self, name: str = "C"):
        """
        @brief Erstellt eine leere Kategorie.
        @param name Bezeichner der Kategorie
        @lastModified 2026-03-10
        """
        # Name der Kategorie
        self.name = name
        # Menge aller Objekte
        self.objects: list[Object] = []
        # Menge aller Morphismen (als Liste)
        self.morphisms: list[Morphism] = []
        # Identitätsmorphismen je Objekt: {obj_name: Morphism}
        self._identities: dict[str, Morphism] = {}

    def add_object(self, obj: Object) -> None:
        """
        @brief Fügt ein Objekt zur Kategorie hinzu.
        @description
            Fügt das Objekt zur Kategorie hinzu und erstellt automatisch
            den Identitätsmorphismus id_obj: obj → obj.
        @param obj Das hinzuzufügende Objekt
        @lastModified 2026-03-10
        """
        # Prüfe, ob Objekt bereits vorhanden
        if obj not in self.objects:
            self.objects.append(obj)
            # Automatisch Identitätsmorphismus erstellen
            identity = Morphism(
                source=obj,
                target=obj,
                name=f"id_{obj.name}",
                func=lambda x: x  # Identitätsfunktion
            )
            self._identities[obj.name] = identity
            self.morphisms.append(identity)

    def add_morphism(self, morphism: Morphism) -> None:
        """
        @brief Fügt einen Morphismus zur Kategorie hinzu.
        @description
            Prüft, dass Quell- und Zielobjekt in der Kategorie existieren,
            dann fügt den Morphismus zur Morphismenmenge hinzu.
        @param morphism Der hinzuzufügende Morphismus
        @raises ValueError falls Quelle oder Ziel nicht in der Kategorie
        @lastModified 2026-03-10
        """
        # Prüfe, ob Quell- und Zielobjekt vorhanden sind
        if morphism.source not in self.objects:
            raise ValueError(f"Quellobjekt {morphism.source.name!r} nicht in Kategorie {self.name!r}")
        if morphism.target not in self.objects:
            raise ValueError(f"Zielobjekt {morphism.target.name!r} nicht in Kategorie {self.name!r}")

        # Morphismus zur Liste hinzufügen (falls nicht bereits vorhanden)
        if morphism not in self.morphisms:
            self.morphisms.append(morphism)

    def compose(self, f: Morphism, g: Morphism) -> Morphism:
        """
        @brief Komposition g ∘ f in der Kategorie.
        @description
            Berechnet die Komposition g ∘ f: A → C für f: A → B und g: B → C.
            Das Ergebnis ist wieder ein Morphismus in der Kategorie.

            Mathematisch: (g ∘ f)(x) = g(f(x))
        @param f Morphismus f: A → B (wird zuerst angewendet)
        @param g Morphismus g: B → C (wird danach angewendet)
        @return Kompositionsmorphismus g ∘ f: A → C
        @raises ValueError falls die Komposition nicht definiert ist
        @lastModified 2026-03-10
        """
        return g.compose(f)

    def identity(self, obj: Object) -> Morphism:
        """
        @brief Gibt den Identitätsmorphismus für ein Objekt zurück.
        @description
            Für jedes Objekt A in einer Kategorie existiert ein eindeutiger
            Identitätsmorphismus id_A: A → A mit:
            - id_A ∘ f = f  (für alle f: X → A)
            - g ∘ id_A = g  (für alle g: A → Y)
        @param obj Das Objekt, für das die Identität gesucht wird
        @return Identitätsmorphismus id_obj
        @raises KeyError falls das Objekt nicht in der Kategorie ist
        @lastModified 2026-03-10
        """
        if obj.name not in self._identities:
            raise KeyError(f"Objekt {obj.name!r} nicht in Kategorie {self.name!r}")
        return self._identities[obj.name]

    def is_functor_valid(self, F: 'Functor') -> bool:
        """
        @brief Prüft, ob F ein gültiger Funktor aus dieser Kategorie ist.
        @description
            Ein Funktor F: C → D ist gültig, wenn:
            1. F(id_A) = id_{F(A)} für alle A ∈ Ob(C)
            2. F(g ∘ f) = F(g) ∘ F(f) für alle komponierbaren f, g

            Diese Methode prüft die Erhaltung von Identitäten und Komposition.
        @param F Der zu prüfende Funktor
        @return True falls F alle Funktoreigenschaften erfüllt
        @lastModified 2026-03-10
        """
        # Prüfe, ob F die korrekte Quellkategorie hat
        if F.source_cat is not self:
            return False

        # Prüfe Erhaltung von Identitäten
        if not F.preserves_identity():
            return False

        # Prüfe Erhaltung von Kompositionen
        if not F.preserves_composition():
            return False

        return True

    def __repr__(self) -> str:
        """Lesbare Darstellung der Kategorie."""
        return f"Category({self.name!r}, {len(self.objects)} Objekte, {len(self.morphisms)} Morphismen)"


# =============================================================================
# KLASSE: Functor
# =============================================================================

class Functor:
    """
    @brief Kovarianter Funktor zwischen zwei Kategorien.
    @description
        Ein kovarianter Funktor F: C → D besteht aus:
        1. Einer Abbildung auf Objekten: F: Ob(C) → Ob(D)
        2. Einer Abbildung auf Morphismen: F: Mor_C(A,B) → Mor_D(F(A),F(B))

        Funktoreigenschaften:
        - Erhaltung der Identität: F(id_A) = id_{F(A)}
        - Erhaltung der Komposition: F(g ∘ f) = F(g) ∘ F(f)

        Beispiele: Vergissfunktor, Hom-Funktor, freier Funktor
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        name: str,
        source_cat: Category,
        target_cat: Category,
        object_map: dict,
        morphism_map: dict
    ):
        """
        @brief Erstellt einen Funktor.
        @param name         Bezeichner des Funktors
        @param source_cat   Quellkategorie C
        @param target_cat   Zielkategorie D
        @param object_map   Abbildung {Object → Object} auf Objekten
        @param morphism_map Abbildung {Morphism → Morphism} auf Morphismen
        @lastModified 2026-03-10
        """
        # Name des Funktors
        self.name = name
        # Quell- und Zielkategorie
        self.source_cat = source_cat
        self.target_cat = target_cat
        # Abbildung auf Objekten: {Object: Object}
        self.object_map = object_map
        # Abbildung auf Morphismen: {Morphism: Morphism}
        self.morphism_map = morphism_map

    def apply_object(self, obj: Object) -> Object:
        """
        @brief Wendet den Funktor auf ein Objekt an.
        @param obj Objekt aus der Quellkategorie
        @return Bild-Objekt in der Zielkategorie
        @raises KeyError falls obj nicht im object_map
        @lastModified 2026-03-10
        """
        if obj not in self.object_map:
            raise KeyError(f"Objekt {obj.name!r} nicht im Objekt-Map des Funktors {self.name!r}")
        return self.object_map[obj]

    def apply_morphism(self, morphism: Morphism) -> Morphism:
        """
        @brief Wendet den Funktor auf einen Morphismus an.
        @param morphism Morphismus aus der Quellkategorie
        @return Bild-Morphismus in der Zielkategorie
        @raises KeyError falls morphism nicht im morphism_map
        @lastModified 2026-03-10
        """
        if morphism not in self.morphism_map:
            raise KeyError(f"Morphismus {morphism.name!r} nicht im Morphismus-Map des Funktors {self.name!r}")
        return self.morphism_map[morphism]

    def preserves_identity(self) -> bool:
        """
        @brief Prüft, ob der Funktor Identitäten erhält.
        @description
            Prüft für alle Objekte A ∈ Ob(C):
            F(id_A) = id_{F(A)}

            Dies ist eine notwendige Eigenschaft eines Funktors.
        @return True falls F(id_A) = id_{F(A)} für alle A
        @lastModified 2026-03-10
        """
        # Prüfe für jedes Objekt in der Quellkategorie
        for obj in self.source_cat.objects:
            try:
                # Hole Identitätsmorphismus in der Quellkategorie
                id_obj = self.source_cat.identity(obj)
                # Wende Funktor auf Identität an
                f_id = self.morphism_map.get(id_obj)
                if f_id is None:
                    # Identitätsmorphismus nicht im Map → prüfe mit apply
                    continue

                # Hole Identität des Bild-Objekts in der Zielkategorie
                img_obj = self.object_map.get(obj)
                if img_obj is None:
                    return False
                id_img = self.target_cat.identity(img_obj)

                # Vergleiche Namen (strukturelle Gleichheit)
                if f_id.name != id_img.name:
                    return False
            except (KeyError, Exception):
                # Objekt oder Morphismus nicht gefunden
                return False
        return True

    def preserves_composition(self) -> bool:
        """
        @brief Prüft, ob der Funktor Kompositionen erhält.
        @description
            Prüft für alle komponierbaren Paare f: A → B, g: B → C:
            F(g ∘ f) = F(g) ∘ F(f)

            Dies ist die zweite wesentliche Eigenschaft eines Funktors.
        @return True falls F alle Kompositionen erhält
        @lastModified 2026-03-10
        """
        # Iteriere über alle Paare von Morphismen in der Quellkategorie
        for f in self.source_cat.morphisms:
            for g in self.source_cat.morphisms:
                # Prüfe ob f und g komponierbar sind: f: A→B, g: B→C
                if f.target != g.source:
                    continue

                # Berechne g ∘ f in der Quellkategorie
                try:
                    gof = self.source_cat.compose(f, g)
                except ValueError:
                    continue

                # F(g ∘ f) – Funktor auf Komposition anwenden
                f_gof = self.morphism_map.get(gof)

                # F(g) und F(f) – Funktor einzeln anwenden
                f_g = self.morphism_map.get(g)
                f_f = self.morphism_map.get(f)

                if f_gof is None or f_g is None or f_f is None:
                    # Nicht im Map → überspringe (kein Fehler, nur unvollständig)
                    continue

                # F(g) ∘ F(f) in der Zielkategorie berechnen
                try:
                    fg_ff = self.target_cat.compose(f_f, f_g)
                except ValueError:
                    return False

                # Vergleiche: F(g ∘ f) = F(g) ∘ F(f)?
                if f_gof.name != fg_ff.name:
                    return False

        return True

    def __repr__(self) -> str:
        """Lesbare Darstellung des Funktors."""
        return (f"Functor({self.name!r}: {self.source_cat.name} → {self.target_cat.name})")


# =============================================================================
# KLASSE: NaturalTransformation
# =============================================================================

class NaturalTransformation:
    """
    @brief Natürliche Transformation zwischen zwei Funktoren.
    @description
        Eine natürliche Transformation α: F ⇒ G zwischen Funktoren
        F, G: C → D besteht aus einer Familie von Morphismen:
        α_A: F(A) → G(A)  für jedes A ∈ Ob(C)

        Diese müssen die Natürlichkeitsbedingung erfüllen:
        Für jeden Morphismus f: A → B gilt:
        α_B ∘ F(f) = G(f) ∘ α_A

        Das heißt: die Naturquadrate kommutieren:
        F(A) --F(f)-→ F(B)
         |               |
        α_A             α_B
         ↓               ↓
        G(A) --G(f)-→ G(B)
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        name: str,
        source_functor: Functor,
        target_functor: Functor,
        components: dict
    ):
        """
        @brief Erstellt eine natürliche Transformation.
        @param name            Bezeichner der Transformation (z.B. "α", "η")
        @param source_functor  Quell-Funktor F: C → D
        @param target_functor  Ziel-Funktor G: C → D
        @param components      {Object → Morphism} (Komponenten α_A: F(A) → G(A))
        @lastModified 2026-03-10
        """
        # Name der natürlichen Transformation
        self.name = name
        # Quell- und Ziel-Funktor
        self.source_functor = source_functor
        self.target_functor = target_functor
        # Komponenten: {Objekt → Morphismus}
        self.components = components

    def naturality_check(self) -> bool:
        """
        @brief Prüft die Natürlichkeitsbedingung für alle Morphismen.
        @description
            Für jeden Morphismus f: A → B in der gemeinsamen Quellkategorie
            muss das Naturquadrat kommutieren:
            α_B ∘ F(f) = G(f) ∘ α_A

            Diese Methode prüft alle nicht-trivialen Morphismen in der
            Quellkategorie und vergleicht die Komponentennamen.
        @return True falls alle Naturquadrate kommutieren
        @lastModified 2026-03-10
        """
        # Quellkategorie (gemeinsam für beide Funktoren)
        source_cat = self.source_functor.source_cat
        target_cat = self.source_functor.target_cat

        # Prüfe für jeden Morphismus f: A → B in der Quellkategorie
        for f in source_cat.morphisms:
            A = f.source
            B = f.target

            # Hole Komponenten α_A und α_B
            alpha_A = self.components.get(A)
            alpha_B = self.components.get(B)

            if alpha_A is None or alpha_B is None:
                # Komponenten fehlen → nicht vollständig definiert
                continue

            # Hole F(f) und G(f)
            Ff = self.source_functor.morphism_map.get(f)
            Gf = self.target_functor.morphism_map.get(f)

            if Ff is None or Gf is None:
                # Morphismen nicht im Map → überspringe
                continue

            # Berechne α_B ∘ F(f)
            try:
                left = target_cat.compose(Ff, alpha_B)
            except ValueError:
                return False

            # Berechne G(f) ∘ α_A
            try:
                right = target_cat.compose(alpha_A, Gf)
            except ValueError:
                return False

            # Vergleiche: α_B ∘ F(f) = G(f) ∘ α_A?
            if left.name != right.name:
                return False

        return True

    def __repr__(self) -> str:
        """Lesbare Darstellung der natürlichen Transformation."""
        return (f"NaturalTransformation({self.name!r}: "
                f"{self.source_functor.name} ⇒ {self.target_functor.name})")


# =============================================================================
# BEISPIELKATEGORIEN
# =============================================================================

def set_category(sets: list) -> Category:
    """
    @brief Erstellt die Set-Kategorie für gegebene Mengen.
    @description
        Die Set-Kategorie hat Mengen als Objekte und Funktionen zwischen
        Mengen als Morphismen. Für je zwei Mengen A, B wird eine generische
        Funktion (als Lambda) als Morphismus hinzugefügt.

        Mathematisch: Set ist die Kategorie der Mengen und Funktionen.
        Ob(Set) = alle Mengen, Mor(A, B) = alle Funktionen f: A → B
    @param sets Liste von Mengen (Python set oder frozenset)
    @return Category-Objekt der Set-Kategorie
    @lastModified 2026-03-10
    """
    # Neue Kategorie "Set" erstellen
    cat = Category("Set")

    # Objekte aus den Mengen erstellen
    objs = []
    for i, s in enumerate(sets):
        name = f"S{i}"
        obj = Object(name, data=s)
        cat.add_object(obj)
        objs.append(obj)

    # Morphismen: Für jedes geordnete Paar eine generische Abbildung
    for i, src_obj in enumerate(objs):
        for j, tgt_obj in enumerate(objs):
            if i != j:  # Identitäten bereits automatisch hinzugefügt
                # Erstelle Morphismus (repräsentiert Existenz einer Funktion)
                morph = Morphism(
                    source=src_obj,
                    target=tgt_obj,
                    name=f"f_{src_obj.name}_{tgt_obj.name}",
                    func=None  # Konkrete Funktion nicht spezifiziert
                )
                cat.add_morphism(morph)

    return cat


def group_category(groups: list) -> Category:
    """
    @brief Erstellt Grp-Kategorie aus einer Liste von Gruppen-Objekten.
    @description
        Die Kategorie Grp hat Gruppen als Objekte und Gruppenhomomorphismen
        als Morphismen.

        Mathematisch: Grp ist die Kategorie der Gruppen.
        Ob(Grp) = alle Gruppen, Mor(G, H) = alle Gruppenhomomorphismen φ: G → H
    @param groups Liste von Gruppen-Bezeichnern oder Gruppenobjekten
    @return Category-Objekt der Grp-Kategorie
    @lastModified 2026-03-10
    """
    # Neue Kategorie "Grp" erstellen
    cat = Category("Grp")

    # Objekte aus den Gruppen erstellen
    objs = []
    for i, g in enumerate(groups):
        # Gruppenname aus dem Objekt lesen oder generieren
        gname = getattr(g, 'name', None) or f"G{i}"
        obj = Object(gname, data=g)
        cat.add_object(obj)
        objs.append(obj)

    # Morphismen: Für jedes Paar einen Homomorphismus-Platzhalter
    for src_obj in objs:
        for tgt_obj in objs:
            if src_obj != tgt_obj:
                morph = Morphism(
                    source=src_obj,
                    target=tgt_obj,
                    name=f"hom_{src_obj.name}_{tgt_obj.name}",
                    func=None
                )
                cat.add_morphism(morph)

    return cat


def opposite_category(cat: Category) -> Category:
    """
    @brief Erstellt die duale (entgegengesetzte) Kategorie C^op.
    @description
        Die duale Kategorie C^op einer Kategorie C entsteht, indem alle
        Morphismenpfeile umgekehrt werden:
        - Ob(C^op) = Ob(C)
        - Mor_{C^op}(A, B) = Mor_C(B, A)

        Mathematisch: Dualitätsprinzip – jede Aussage in C entspricht
        einer Aussage in C^op mit umgekehrten Pfeilen.
    @param cat Die ursprüngliche Kategorie C
    @return Neue Kategorie C^op
    @lastModified 2026-03-10
    """
    # Neue duale Kategorie erstellen
    op_cat = Category(f"{cat.name}^op")

    # Objekte bleiben gleich
    for obj in cat.objects:
        op_cat.add_object(obj)

    # Morphismen werden umgekehrt (außer Identitäten, die bleiben gleich)
    for morph in cat.morphisms:
        # Identitätsmorphismen bleiben Identitäten
        if morph.source == morph.target and morph.name.startswith("id_"):
            continue  # Schon durch add_object hinzugefügt

        # Morphismus f: A → B wird zu f^op: B → A
        op_morph = Morphism(
            source=morph.target,   # Quelle und Ziel tauschen
            target=morph.source,
            name=f"{morph.name}^op",
            func=None  # Umkehrfunktion nicht immer definierbar
        )
        op_cat.add_morphism(op_morph)

    return op_cat


def product_category(cat1: Category, cat2: Category) -> Category:
    """
    @brief Erstellt die Produktkategorie C₁ × C₂.
    @description
        Die Produktkategorie C₁ × C₂ hat:
        - Objekte: Paare (A, B) mit A ∈ Ob(C₁) und B ∈ Ob(C₂)
        - Morphismen: Paare (f, g): (A,B) → (A',B') mit f: A→A' und g: B→B'
        - Komposition: komponentenweise (f'∘f, g'∘g)
        - Identitäten: (id_A, id_B)
    @param cat1 Erste Kategorie C₁
    @param cat2 Zweite Kategorie C₂
    @return Produktkategorie C₁ × C₂
    @lastModified 2026-03-10
    """
    # Neue Produktkategorie erstellen
    prod_cat = Category(f"{cat1.name}×{cat2.name}")

    # Objekte: Alle Paare (A, B)
    obj_pairs: dict[tuple, Object] = {}
    for obj1 in cat1.objects:
        for obj2 in cat2.objects:
            # Erstelle Paar-Objekt
            pair_name = f"({obj1.name},{obj2.name})"
            pair_obj = Object(pair_name, data=(obj1, obj2))
            prod_cat.add_object(pair_obj)
            obj_pairs[(obj1.name, obj2.name)] = pair_obj

    # Morphismen: Alle Paare (f, g) von Morphismen
    for morph1 in cat1.morphisms:
        for morph2 in cat2.morphisms:
            # Paar-Morphismus (f, g): (A,B) → (A',B')
            src_name = (morph1.source.name, morph2.source.name)
            tgt_name = (morph1.target.name, morph2.target.name)

            src_pair = obj_pairs.get(src_name)
            tgt_pair = obj_pairs.get(tgt_name)

            if src_pair is None or tgt_pair is None:
                continue

            # Überspringe triviale Identitätspaare (bereits automatisch erstellt)
            if src_pair == tgt_pair and morph1.name.startswith("id_") and morph2.name.startswith("id_"):
                continue

            pair_morph = Morphism(
                source=src_pair,
                target=tgt_pair,
                name=f"({morph1.name},{morph2.name})",
                func=None
            )
            prod_cat.add_morphism(pair_morph)

    return prod_cat


def is_isomorphism(f: Morphism, cat: Category) -> bool:
    """
    @brief Prüft ob f ein Isomorphismus in der Kategorie ist.
    @description
        Ein Morphismus f: A → B ist ein Isomorphismus, falls es einen
        Morphismus g: B → A gibt (das Inverse), sodass:
        g ∘ f = id_A  und  f ∘ g = id_B

        In Set: Bijektive Funktionen.
        In Grp: Bijektive Gruppenhomomorphismen.
    @param f   Der zu prüfende Morphismus
    @param cat Die Kategorie, in der geprüft wird
    @return True falls f ein Isomorphismus ist
    @lastModified 2026-03-10
    """
    A = f.source
    B = f.target

    # Suche Inverse: g: B → A mit g ∘ f = id_A und f ∘ g = id_B
    for g in cat.morphisms:
        if g.source != B or g.target != A:
            continue

        # Prüfe g ∘ f = id_A
        try:
            gof = cat.compose(f, g)
            id_A = cat.identity(A)
            if gof.name != id_A.name:
                continue
        except (ValueError, KeyError):
            continue

        # Prüfe f ∘ g = id_B
        try:
            fog = cat.compose(g, f)
            id_B = cat.identity(B)
            if fog.name != id_B.name:
                continue
        except (ValueError, KeyError):
            continue

        return True  # Inverse g gefunden

    return False


def is_epimorphism(f: Morphism, cat: Category) -> bool:
    """
    @brief Prüft ob f ein Epimorphismus ist (rechts-kürzbar).
    @description
        Ein Morphismus f: A → B ist ein Epimorphismus, falls für alle
        Morphismen g, h: B → C gilt:
        g ∘ f = h ∘ f  ⟹  g = h

        In Set: Surjektive Funktionen.
        Diese Methode prüft die Kürzungseigenschaft für alle vorhandenen
        Morphismenpaare in der Kategorie.
    @param f   Der zu prüfende Morphismus
    @param cat Die Kategorie
    @return True falls f ein Epimorphismus ist
    @lastModified 2026-03-10
    """
    B = f.target

    # Finde alle Paare (g, h): B → C
    for g in cat.morphisms:
        if g.source != B:
            continue
        for h in cat.morphisms:
            if h.source != B or h.target != g.target:
                continue
            if g == h:
                continue  # Triviale Gleichheit

            # Prüfe ob g ∘ f = h ∘ f
            try:
                gof = cat.compose(f, g)
                hof = cat.compose(f, h)
            except ValueError:
                continue

            # Wenn g ∘ f = h ∘ f aber g ≠ h → kein Epimorphismus
            if gof.name == hof.name and g.name != h.name:
                return False

    return True


def is_monomorphism(f: Morphism, cat: Category) -> bool:
    """
    @brief Prüft ob f ein Monomorphismus ist (links-kürzbar).
    @description
        Ein Morphismus f: A → B ist ein Monomorphismus, falls für alle
        Morphismen g, h: C → A gilt:
        f ∘ g = f ∘ h  ⟹  g = h

        In Set: Injektive Funktionen.
        Diese Methode prüft die Kürzungseigenschaft für alle vorhandenen
        Morphismenpaare in der Kategorie.
    @param f   Der zu prüfende Morphismus
    @param cat Die Kategorie
    @return True falls f ein Monomorphismus ist
    @lastModified 2026-03-10
    """
    A = f.source

    # Finde alle Paare (g, h): C → A
    for g in cat.morphisms:
        if g.target != A:
            continue
        for h in cat.morphisms:
            if h.target != A or h.source != g.source:
                continue
            if g == h:
                continue  # Triviale Gleichheit

            # Prüfe ob f ∘ g = f ∘ h
            try:
                fog = cat.compose(g, f)
                foh = cat.compose(h, f)
            except ValueError:
                continue

            # Wenn f ∘ g = f ∘ h aber g ≠ h → kein Monomorphismus
            if fog.name == foh.name and g.name != h.name:
                return False

    return True


def hom_set(A: Object, B: Object, cat: Category) -> list:
    """
    @brief Gibt alle Morphismen von A nach B zurück (Hom(A, B)).
    @description
        Das Hom-Set Hom_C(A, B) ist die Menge aller Morphismen von A nach B:
        Hom_C(A, B) = {f ∈ Mor(C) | f: A → B}

        In Set: Hom_Set(A, B) = B^A (Menge aller Funktionen von A nach B).
    @param A   Quellobjekt
    @param B   Zielobjekt
    @param cat Die Kategorie
    @return Liste aller Morphismen von A nach B
    @lastModified 2026-03-10
    """
    # Filtere alle Morphismen von A nach B
    return [m for m in cat.morphisms if m.source == A and m.target == B]


def yoneda_functor(cat: Category, obj: Object) -> Functor:
    """
    @brief Yoneda-Einbettung: h_A = Hom(A, -) als Funktor.
    @description
        Der Yoneda-Funktor h_A: C → Set ist definiert durch:
        - Auf Objekten: h_A(B) = Hom_C(A, B) (als Objekt der Set-Kategorie)
        - Auf Morphismen: h_A(f) = f ∘ - (Nachkomposition mit f)

        Das Yoneda-Lemma besagt:
        Nat(h_A, F) ≅ F(A)  für jeden Funktor F: C → Set

        Die Einbettung A ↦ h_A ist vollständig und treu (Yoneda-Einbettung).
    @param cat Die Quellkategorie C
    @param obj Das Objekt A, für das h_A konstruiert wird
    @return Funktor h_A: C → Set
    @lastModified 2026-03-10
    """
    # Zielkategorie: Set (Mengen von Morphismen als Objekte)
    set_cat = Category(f"Set[{cat.name}]")

    # Erstelle Objekte in Set: Für jedes B ∈ Ob(C) wird Hom(A, B) ein Objekt
    obj_map = {}
    for B in cat.objects:
        # Das Hom-Set Hom(A, B) wird als Menge von Morphismennamen repräsentiert
        hom_AB = hom_set(obj, B, cat)
        hom_set_name = f"Hom({obj.name},{B.name})"
        hom_obj = Object(hom_set_name, data=hom_AB)
        set_cat.add_object(hom_obj)
        obj_map[B] = hom_obj

    # Erstelle Morphismen in Set: Für jeden f: B → C wird h_A(f): Hom(A,B) → Hom(A,C)
    morph_map = {}
    for f in cat.morphisms:
        B = f.source
        C_obj = f.target

        src_hom = obj_map.get(B)
        tgt_hom = obj_map.get(C_obj)

        if src_hom is None or tgt_hom is None:
            continue

        # Nachkomposition: h_A(f)(g) = f ∘ g für g: A → B
        f_morph = Morphism(
            source=src_hom,
            target=tgt_hom,
            name=f"Hom({obj.name},{f.name})",
            func=None  # Repräsentation der Nachkomposition
        )
        set_cat.add_morphism(f_morph)
        morph_map[f] = f_morph

    return Functor(
        name=f"h_{obj.name}",
        source_cat=cat,
        target_cat=set_cat,
        object_map=obj_map,
        morphism_map=morph_map
    )
