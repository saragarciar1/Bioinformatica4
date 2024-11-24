#Subir los datos
datos <- read.table("C:\\Users\\sgarcrot\\Downloads\\datostrabajoR.txt", header=TRUE)

#Analizar los datos
head (datos) #para mostrar las primeras filas
summary(datos) #resumen de cada columna
dim(datos) #dimension de las filas y columnas
 str(datos)#estructura de los datos
 
 #Numero de variables
 num_variables <-ncol(datos) 
 #numero de tratamientos
num_tratamiento <- length(unique(datos$tratamiento)) 
datos$tratamiento <- as.factor(datos$tratamiento)
cat("Numero de variables:", num_variables, "\n")
cat("Numero de tratamientos:", num_tratamiento, "\n")

 #elegir colores para cada condicion
 colores_condiciones <- c("Wildtype"="red", "Sequia"="pink", "ExcesoRiego"="green")

 #hacer el boxplot
 boxplot(datos$wiltdtype, datos$sequia, datos$excesoriego,
 	names = c("wildtype", "sequia", "excesoriego"),
 	col = colores_condiciones,
 	main ="Boxplot de cada condición",
 	ylab ="valores de Medición",
 	xlab="condiciones")
 
 #definir una paleta de colores para los tratamiento
 colores_tratamiento <-rainbow (length(unique(datos$tratamiento)))
 
 #grafico de sequia vs. wildtype
 plot(datos$wildtype, datos$sequia,
 	col=colores_tratamiento[datos$tratamiento],
 	pch=19,
 	xlab="wildtype",
 	ylab="sequia",
 	main="dispersion de sequia vs wildtype")
 	
 #leyenda en la esquina inferior derecha
 legend("bottomright", legend=unique(datos$tratamiento),
 	col=colores_tratamiento, pch=19, title="tratamiento")
 
 #grafico de excesoriego VS wildtype
 plot(datos$wildtype, datos$excesoriego,
 	col=colores_tratamiento[datos$tratamiento],
 	pch=19,
 	xlab="wildtype",
 	ylab="excesoriego",
 	main="dispersion de excesoriego vs wildtype")
 	
 #leyenda en la esquina inferior derecha
 legend ("bottomright", legend=unique(datos$tratamiento),
 	col=colores_tratamiento, pch=19, title="tratamiento")
 
 #histograma para wildtype
 hist(datos$wildtype, main="histograma de wildtype", col=colores_condiciones["wildtype"], xlab="wildtype")
 
 #histograma para sequia
 hist(datos$sequia, main="histograma de sequia", col=colores_condiciones["sequia"], xlab="sequia")
 
 #hisograma para excesoriego
 hist(datos$excesoriego, main="histograma de excesoriego", col=colores_condiciones["excesoriego"], xlab="excesoriego")
 
 #convertir la columna tratamiento a factor
 tratamiento_factor <-factor(datos$tratamiento)
 
 #media por tratamiento
 medias <- aggregate(cbind(wildtype, sequia, excesoriego) ~ tratamiento, data=datos, FUN = mean)
 
 #desviacion estandar por tratamiento
 desviaciones <- aggregate(cbind(wildtype, sequia, excesoriego) ~ tratamiento, data=datos, FUN=sd)
 
 print(medias)
 print(desviaciones)
 
 #contar elementos en cada tratamiento
 table(tratamiento_factor)
 
 #extraer datos para tratamientos 1 y 4
 datos_tratamiento1 <- subset(datos, tratamiento==1)
 datos_tratamiento4 <- subset(datos, tratamiento==4)
 
 #prueba de normalidad para wildtype y sequia en tratamiento
 datos_tratamiento5 <- subset(datos, tratamiento==5)
 shapiro.test(datos_tratamiento1$wildtype)
 shapiro.test(datos_tratamiento1$sequia)
 shapiro.test(datos_tratamiento5$wildtype)
 shapiro.test(datos_tratamiento5$sequia)
 
 #compraracion entre wildtype y sequia para tratamiento 1 y 5 
 t.test(datos_tratamiento1$wildtype, datos_tratamiento1$sequia)
 t.test(datos_tratamiento5$wildtype, datos_tratamiento5$sequia)
 
 #compraracion entre sequia y excesoriego
 t.test(datos_tratamiento1$sequia, datos_tratamiento1$excesoriego)
 t.test(datos_tratamiento5$sequia, datos_tratamiento5$excesoriego)
 
 #separar las condiciones de tratamiento 1
 wildtype_tratamiento1 <- datos_tratamiento1$wildtype
 sequia_tratamiento1 <- datos_tratamiento1$sequia
 excesoriego_tratamiento1 <- datos_tratamiento1$excesoriego
 
 #crear tabla para ANOVA
 valores <- c(wildtype_tratamiento1, sequia_tratamiento1, excesoriego_tratamiento1)
 condicion <- rep(c("wildtype", "sequia", "excesoriego"),
                  each=length(wildtype_tratamiento1))
 datos_anova <- data.frame(
 	condicion=condicion,
 	valor=valores
 )
 #ejecutar el ANOVA
 anova_resultado <- aov(valor ~ condicion, data=datos_anova)
summary(anova_resultado) 
 
