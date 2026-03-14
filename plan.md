# Plan de Trabajo: Física Estadística Avanzada - Teoría Cinética de Sistemas Pequeños

## Fase 1: Parte Teórica

### 1.1 Colectividad Microcanónica
- [ ] Derivar función densidad ρ(p) integrando sobre todos los grados de libertad excepto el momento de una partícula
- [ ] Obtener distribución de velocidades v
- [ ] Obtener distribución de energía ε
- [ ] Tomar límite N→∞ y verificar que se recupera Maxwell-Boltzmann

### 1.2 Colectividad de Dinámica Molecular
- [ ] Escribir la función densidad dada: ρ = δ(P - Σp_i) e^{-βH}
- [ ] Obtener función densidad para momento de una partícula con P=0
- [ ] Comparar con resultado microcanónico

### 1.3 Método de la Transformada de Laplace
- [ ] Aplicar transformada de Laplace respecto a la energía
- [ ] Realizar integrales sobre momentos
- [ ] Aplicar transformada inversa para obtener distribuciones
- [ ] Repetir para colectividad microcanónica y dinámica molecular
- [ ] En DM: aplicar transformada de Fourier respecto al momento total

---

## Fase 2: Parte Práctica - Programa v1 (Microcanónica)

### 2.1 Estructura básica
- [x] Crear programa con arrays vx(N), vy(N), vz(N) para N≤32 partículas
- [x] Inicializar velocidades con igual módulo y dirección aleatoria
- [x] Imponer momento total cero: restar velocidad del centro de masas

### 2.2 Control térmico
- [x] Escalar velocidades para fijar temperatura: v_i' = v_i √(T_deseada/T_actual)

### 2.3 Colisiones partículas-partícula
- [x] Seleccionar pares de partículas aleatoriamente
- [x] Implementar colisión elástica 2D/3D
- [x] Actualizar velocidades de ambas partículas

### 2.4 Condiciones de contorno periódicas
- [x] Reemplazar colisiones con pared por PBC
- [x] Cuando partícula sale por un lado, entra por el opuesto
- [x] Mantener conservación exacta del momento total

### 2.5 Muestreo y análisis
- [x] Después de ~10 colisiones/partícula, muestrear distribución
- [x] Repetir muestreo acumulativamente
- [x] Generar histograma y comparar con Maxwell-Boltzmann teórica

---

## Fase 3: Parte Práctica - Programa v2 (Dinámica Molecular)

### 3.1 Adaptación del programa v1
- [ ] Copiar programa v1
- [ ] Eliminar colisiones con paredes
- [ ] Mantener solo colisiones elásticas entre partículas

### 3.2 Verificación
- [ ] Verificar conservación de momento total
- [ ] Verificar conservación de energía
- [ ] Comparar distribución de velocidades con teoría microcanónica

---

## Fase 4: Parte Práctica - Programa v3 (DSMC)

### 4.1 Método DSMC
- [ ] Buscar paper de Montanero y Santos sobre DSMC
- [ ] Implementar algoritmo DSMC completo
- [ ] Muestrear múltiples sistemas simultáneamente

### 4.2 Comparación
- [ ] Comparar resultados v3 con v1 y v2
- [ ] Comparar con predicciones teóricas de ambas colectividades

---

## Fase 5: Documentación y Entrega

- [ ] Escribir memoria explicando teoría y resultados
- [ ] Incluir gráficos comparativos teoría vs simulación
- [ ] Documentar código

---

## Notas Técnicas

### Entorno recomendado
- Lenguaje: Python (numpy, matplotlib) o C/Fortran
- N partículas: ≤32 (según especificación)
- Colisiones típicas: ~10 por partícula para equilibrio

### Verificaciones importantes
- Momento total = 0 (tras inicialización)
- Energía constante (en microcanónica)
- Distribución converge a Maxwell-Boltzmann en límite N→∞
