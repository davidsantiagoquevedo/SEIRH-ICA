# SEIRHT-ICA

En este repositorio se presentan los archivos asociados al proyecto de modelo determinista tipo SEIR con estados hospitalarios y aislamiento de infecciosos  basado en el proyecto: https://github.com/judmejiabe/SEIIHR/.
# Índice de documentos
* `ModeloConAislamientoICA.pdf`: Documento técnico con descripción del modelo de compartimentos, el sistema de ecuaciones diferenciales asociados al mismo y resultados preliminares para ciudad de Bogotá con distintos escenarios de testeo y cuarentena.
* `SEIRH_Class_ICA.py`: Script de python con solución numérica del modelo determinista.

# Descripción del modelo de compartimentos
![Esquema](/IMG/Esquema.png?raw=true)

* S: Susceptibles
* E: Expuestos
* I<sub>M</sub>: Infectados moderados
* I<sub>C</sub>: Infectados que requerirán algún tipo de hospitalización pero aún no se han aislado en su casa
* I<sub>CA</sub>: Infectados que requerirán y se encuentran aislados en su casa
* I<sub>HR</sub>: Infectados hospitalizados que se recuperan
* I<sub>UR</sub>: Infectados en UCI que se recuperan
* I<sub>HD</sub>: Infectados hospitalizados que fallecen
* I<sub>UD</sub>: Infectados en UCI que fallecen

* I<sub>R</sub>: Infectados que pasan a camilla de recuperación después de estar en UCI
* R: Recuperados
* D: Fallecidos

# Tasas de transición entre compartimentos
![Tasas](/IMG/Tasas.png?raw=true)

# Sistema de ecuaciones diferenciales ordinarias
![Ecuaciones](/IMG/Ec.png?raw=true)

                
