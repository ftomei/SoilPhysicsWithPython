# Soil Water Content, Porosity, and Saturation

The functions `computeporosity()` and `ComputeSaturationWetness()` calculate several important soil physical properties from bulk density and gravimetric water content.

This documentation explains the physical relationships and equations used by these functions and shows how they are implemented in Python.

---

## 1. Input Variables and Constants

The functions use the following quantities:

| Symbol | Variable                  | Description               | Unit  |
| ------ | ------------------------- | ------------------------- | ----- |
| `ρb`   | `bulkdensity`             | Soil bulk density         | kg/m³ |
| `ρp`   | `ParticleDensity`         | Soil particle density     | kg/m³ |
| `ρw`   | `waterDensity`            | Water density             | kg/m³ |
| `w`    | `gravWaterContent`        | Gravimetric water content | kg/kg |
| `θ`    | `WaterContent_Volumetric` | Volumetric water content  | m³/m³ |
| `φ`    | `porosity`                | Total porosity            | —     |
| `e`    | `voidRatio`               | Void ratio                | —     |
| `S`    | `degree_of_saturation`    | Degree of saturation      | —     |
| `φg`   | `GasPorosity`             | Gas-filled porosity       | —     |


The code assumes the following constant values:

```python
waterDensity = 1000
ParticleDensity = 2650
```

Therefore:

**ρw = 1000 kg/m³**

and:

**ρp = 2650 kg/m³**

The particle density is assumed to be constant for these calculations.

---

## 2. Gravimetric and Volumetric Water Content

### Gravimetric Water Content

Gravimetric water content is the ratio of the mass of water to the mass of dry soil:

**w = Mw / Ms**

where:

* **Mw** is the mass of water.
* **Ms** is the mass of dry soil.

The units are usually expressed as kg/kg.

### Volumetric Water Content

Volumetric water content is the volume of water per total volume of soil:

**θ = Vw / Vt**

The relationship between gravimetric and volumetric water content is:

**θ = w × (ρb / ρw)**

where:

* **θ** = volumetric water content.
* **w** = gravimetric water content.
* **ρb** = bulk density.
* **ρw** = water density.

This relationship is implemented in `computeporosity()` as:

```python
WaterContent_Volumetric = gravWaterContent * (
    bulkdensity / waterDensity
)
```

Therefore, the input `gravWaterContent` is converted into volumetric water content using the ratio of bulk density to water density.

---

## 3. Total Porosity

Total porosity is the fraction of the total soil volume that is occupied by pore space.

It can be calculated from bulk density and particle density as:

**φ = 1 − (ρb / ρp)**

where:

* **φ** = total porosity.
* **ρb** = bulk density.
* **ρp** = particle density.

The equation is implemented as:

```python
porosity = 1 - (bulkdensity / ParticleDensity)
```

For example, if:

**ρb = 1300 kg/m³**

and:

**ρp = 2650 kg/m³**

then:

**φ = 1 − (1300 / 2650)**

**φ ≈ 0.509**

Therefore, approximately **50.9%** of the soil volume consists of pore space.

---

## 4. Void Ratio

Void ratio is the ratio of the volume of pores to the volume of solids:

**e = Vv / Vs**

Porosity and void ratio are related by:

**e = φ / (1 − φ)**

The Python implementation is:

```python
voidRatio = porosity / (1 - porosity)
```

Porosity and void ratio are different quantities:

* **Porosity** is the fraction of the total soil volume occupied by pores.
* **Void ratio** is the ratio of pore volume to solid volume.

---

## 5. Volumetric Water Content at Saturation

When a soil is completely saturated, all of the pore space is filled with water.

Therefore:

**Vw = Vv**

As a result, the volumetric water content at saturation is equal to the total porosity:

**θsat = φ**

This relationship is used to calculate the gravimetric water content at saturation.

---

## 6. Gravimetric Water Content at Saturation

The general relationship between volumetric and gravimetric water content is:

**θ = w × (ρb / ρw)**

At saturation:

**θsat = ws × (ρb / ρw)**

Since:

**θsat = φ**

we obtain:

**φ = ws × (ρb / ρw)**

Solving for the gravimetric water content at saturation gives:

**ws = φ / (ρb / ρw)**

The total porosity is:

**φ = 1 − (ρb / ρp)**

Substituting the porosity equation gives:

**ws = [1 − (ρb / ρp)] / (ρb / ρw)**

This is the equation used by the `ComputeSaturationWetness()` function.

The corresponding Python implementation is:

```python
porosity = 1 - (bulkdensity / ParticleDensity)

gravWaterContent_at_saturation = (
    porosity / (bulkdensity / waterDensity)
)
```

The function first calculates total porosity and then converts the saturated volumetric water content into its corresponding gravimetric water content.

---

## 7. Degree of Saturation

The degree of saturation represents the fraction of the pore space that is filled with water.

It is defined as:

**S = θ / φ**

where:

* **S** = degree of saturation.
* **θ** = volumetric water content.
* **φ** = total porosity.

The Python implementation is:

```python
degree_of_saturation = WaterContent_Volumetric / porosity
```

At complete saturation:

**θ = φ**

Therefore:

**S = φ / φ = 1**

A value of **S = 1** means that all pore space is filled with water.

---

## 8. Gas Porosity

Gas porosity represents the fraction of the total soil volume occupied by air or gas.

The total pore space consists of water-filled and gas-filled pores:

**φ = θ + φg**

where **φg** is gas porosity.

Therefore:

**φg = φ − θ**

The Python implementation is:

```python
GasPorosity = porosity - WaterContent_Volumetric
```

At saturation:

**θ = φ**

Therefore:

**φg = 0**

This means that there is no air-filled pore space when the soil is completely saturated.

---

## 9. Relationship Between the Calculated Properties

The main relationships used by the functions can be summarized as follows.

### Total Porosity

**φ = 1 − (ρb / ρp)**

### Void Ratio

**e = φ / (1 − φ)**

### Volumetric Water Content

**θ = w × (ρb / ρw)**

### Gas Porosity

**φg = φ − θ**

### Degree of Saturation

**S = θ / φ**

### Volumetric Water Content at Saturation

**θsat = φ**

### Gravimetric Water Content at Saturation

**ws = φ / (ρb / ρw)**

or, after substituting the porosity equation:

**ws = [1 − (ρb / ρp)] / (ρb / ρw)**

---

## 10. Worked Example

Consider a soil with:

**ρb = 1300 kg/m³**

and:

**w = 0.20 kg/kg**

using:

**ρp = 2650 kg/m³**

and:

**ρw = 1000 kg/m³**

### Step 1: Calculate Total Porosity

**φ = 1 − (1300 / 2650)**

**φ ≈ 0.509**

### Step 2: Calculate Volumetric Water Content

**θ = 0.20 × (1300 / 1000)**

**θ = 0.260**

### Step 3: Calculate Gas Porosity

**φg = 0.509 − 0.260**

**φg ≈ 0.249**

### Step 4: Calculate Degree of Saturation

**S = 0.260 / 0.509**

**S ≈ 0.511**

Therefore, the soil is approximately **51.1% saturated**.

### Step 5: Calculate Gravimetric Water Content at Saturation

At saturation:

**θsat = φ ≈ 0.509**

Then:

**ws = 0.509 / (1300 / 1000)**

**ws ≈ 0.391 kg/kg**

Therefore, the gravimetric water content required for saturation is approximately **0.391 kg/kg**.

Since the input water content was:

**w = 0.20 kg/kg**

and:

**0.20 < 0.391**

the input water content is consistent with a soil that is not fully saturated.

---

## 11. Relationship to the Python Functions

The `computeporosity()` function performs the following calculations:

1. Total porosity.
2. Void ratio.
3. Volumetric water content.
4. Gas porosity.
5. Degree of saturation.

The `ComputeSaturationWetness()` function calculates the gravimetric water content corresponding to complete saturation:

```python
def ComputeSaturationWetness(bulkdensity):
    waterDensity = 1000
    ParticleDensity = 2650

    porosity = 1 - (bulkdensity / ParticleDensity)

    gravWaterContent_at_saturation = (
        porosity / (bulkdensity / waterDensity)
    )

    return gravWaterContent_at_saturation
```

The returned value can then be compared with the supplied gravimetric water content.

In `main()`:

```python
if (setMassWetness > gravWaterContent >= 0):
    computeporosity(bulkDensity, gravWaterContent)
else:
    print("Invalid input")
```

This check ensures that the supplied gravimetric water content is non-negative and does not exceed the calculated gravimetric water content at saturation.

---

## 12. Summary

The calculations in these functions are based on a small set of fundamental soil-water relationships.

### Volumetric Water Content

**θ = w × (ρb / ρw)**

### Total Porosity

**φ = 1 − (ρb / ρp)**

### Volumetric Water Content at Saturation

**θsat = φ**

### Gravimetric Water Content at Saturation

**ws = [1 − (ρb / ρp)] / (ρb / ρw)**

Together, these relationships connect the physical concepts of **bulk density, particle density, porosity, water content, gas-filled pore space, and degree of saturation** with their implementation in Python.
