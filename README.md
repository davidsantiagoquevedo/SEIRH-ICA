# SEIRH-ICA

En este repositorio se presentan los archivos asociados al proyecto de modelo determinista tipo SEIR con estados hospitalarios y aislamiento de infecciosos  basado en el proyecto: https://github.com/judmejiabe/SEIIHR/.
# Índice de documentos
* `ModeloConAislamientoICA.pdf`: Documento técnico con descripción del modelo de compartimentos, el sistema de ecuaciones diferenciales asociados al mismo y resultados preliminares para ciudad de Bogotá comparados con el escenario con A=60% y R0=3.0 del modelo https://github.com/judmejiabe/SEIIHR/.
* `SEIRH_Class_ICA.py`: Script de python con solución numérica del modelo determinista.
# Estructura del Código Principal `SEIRH_Class_ICA.py`
La estructura del código principal se desarrolla en torno a la clase `SEIRH()`, la cual está compuesta de las siguientes funciones:
* `var_trans(self,beta_0,beta_1,beta_H,r,A,Movilidad)`: Permite definir las variables de transmisión fijas del modelo.
* `beta(self,t)`: Calcula el valor de Beta dependiendo de las condiciones de  movilidad.
* `var_t_estadia(self, omega, gamma_M, sigma_C, sigma_CA, gamma_HR, nu, gamma_R, sigma_HD, sigma_UD)`: Fija las variables de transición entre compartimentos, inversas al tiempo de estadía en estos.
* `var_H_UCI(self, delta_M, delta_HR, delta_HD, delta_UR)`: Define las probabilidades de transición a cada estado hospitalario.
* `var_ini(self, N0, E0=0, IM0=0, IC0=0, ICA0=0,IHR0=0, IUR0=0, IHD0=0, IUD0=0, IR0=0, R0=0, D0=0 )`: Fija las condiciones iniciales del modelo.
* `ODES(self,y,t)`: Define el sistema de ecuaciones diferenciales asociado al modelo.
* `solve(self,t0,tf)`: Utiliza la función  `odeint()` de python para resolver numéricamente el sistema de ecuaciones diferenciales.

# Ejecutar simulaciones
Para ejecutar las simulaciones se utiliza la función `solve()`. Es importante tener en cuenta que el tiempo en que se comienza a correr depende del escenario que se busque reproducir, en el caso en que se ejecuten los escenarios que comienzan a partir del 12 de abril, se debe asignar `t0=42`, dado que este es el día 42 desde el inicio de la pandemia, el cual se cuenta con respecto al primer dato de la base de movilidad (marzo 1 de 2020).

Cabe resaltar que para ejecutar las simulaciones es necesario definir los atributos del modelo con las funciones 1 a 5 de la lista anterior.

# Bases de Datos
Se encuentran guardadas en el directorio BBDD, se listan a continuación:
* `A.xlsx`: Índice de movilidad desde marzo 1 hasta junio 1
* `RES_SecSalud_Rev`: Resultados de los escenarios de la secretaría de salud sin considerar aislamiento de individuos ICA.

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

                
