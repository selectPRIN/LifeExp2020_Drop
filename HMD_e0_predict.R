rm(list=ls())
stmf <- read.csv("https://www.mortality.org/Public/STMF/Outputs/stmf.csv",header=TRUE,sep=",",skip=2)
levels(stmf$CountryCode)[1] <- "AUS"
#Function to read HMD life table
hmd.lt <- function (country, username, password, label = country) 
{
    path <- paste("https://www.mortality.org/hmd/", country, 
        "/STATS/", "bltper_1x1.txt", sep = "")
    userpwd <- paste(username, ":", password, sep = "")
    txt <- RCurl::getURL(path, userpwd = userpwd)
    con <- textConnection(txt)
    life.table <- try(utils::read.table(con, skip = 2, header = TRUE, 
        na.strings = "."), TRUE)
    close(con)
    if (class(life.table) == "try-error") 
        stop("Connection error at www.mortality.org. Please check username, password and country label.")
    return(life.table)
}
nax <- function(lt,age){  # determine nax from HMD 1x1 LT for age interval [x1,xn]
    # Credits to Carl Boe
  j1=age[-length(age)]+1; jn=age[-1]+1;
  nLx =( lt$T[j1] - lt$T[jn]);
  ndx=lt$l[j1]-lt$l[jn] ;
  nax = (nLx - (jn-j1)* lt$l[jn])/ ndx;
  return(nax)
  }

get.lt <- function(mx,age,sex="F",ax){#Simple function to get life table from rates, for simplicity ax are imported from HMD tables
    ageint <- age[-1]-age[-length(age)]
 #   ax <- rep(NA,length(age))
 #   ax[1:length(age)-1] <- ageint/2
    ax <- c(ax,ifelse(mx[length(age)]>0,1/mx[length(age)],0))
    qx <- rep(NA,length(age))
    qx[1:length(age)-1] <- ageint*mx[1:length(age)-1]/(1+(ageint-ax[1:length(age)-1])*mx[1:length(age)-1])
    qx[length(age)] <- 1
    lx <- rep(NA,length(age))
    lx <- round(c(100000,100000*cumprod(1-qx[-length(age)])))
    lx[lx<0|is.na(lx)] <- 0
    dx <- c(-diff(lx),lx[length(lx)])
    Lx <- (c(ageint,110-age[length(age)]))*lx-(c(ageint,110-age[length(age)])-ax)*dx
    Lx[length(age)] <- sum(c(0,lx[length(age)]*ax[length(age)]),na.rm=TRUE)
    Tx <- rev(cumsum(rev(Lx)))
    ex <- ifelse(lx > 0, Tx/lx , NA)
    
    return(data.frame(x=age,nax=ax,nmx=mx,nqx=qx,lx=lx,
                     ndx=dx,nLx=Lx,Tx=Tx,ex=ex))
}

age <- c(0,15,65,75,85)

mymatrix <- data.frame("Country"=NA,"Year"=NA,"e0.smtf"=NA,"e0.hmd"=NA,"ratio014"=NA,"ratio1564"=NA,"ratio6574"=NA,"ratio7584"=NA,"ratio85p"=NA)

for (i in 1:length(names(table(stmf$CountryCode)))){

    
    cnt <- names(table(stmf$CountryCode))[i]

    #get input data, to include deaths with missing dates 
    input <- read.csv(paste("STMFinput/",cnt,"stmf.csv",sep=""),header=TRUE,sep=",",na.strings=".")
        
        #GET life table from HMD
    cnt.lt <- hmd.lt(country=cnt,user="USERNAME",password="PASSWORD")

    # Compute life expectancy up to 2019 (using both stmf and hmd data
    for(y in as.numeric(names((table(stmf$Year[stmf$CountryCode == cnt]))))[1]:2019){
        D.y <- sum(stmf[which(stmf$CountryCode==cnt&stmf$Sex=="b"&stmf$Year==y),10],na.rm=TRUE)
        D.nodate <- input$Deaths[which(input$Sex=="b"&input$Year==y&input$Week=="UNK"&input$Age=="TOT")]
        D.nodate <- ifelse(length(D.nodate)==0,0,D.nodate)

        #RATES. Note I adjust rates by including deaths with missing dates. I'm assuming that deaths with mssing dates are uniformly distribued across age classes,
        #There are some countries (e.g. Sweden) where the number of deaths with missing dates is relatively high
        mx.y <- ((D.nodate+D.y)/D.y)*colSums(stmf[which(stmf$CountryCode==cnt&stmf$Sex=="b"&stmf$Year==y),11:15])/52

        m.hmd <- c(sum(cnt.lt$dx[cnt.lt$Year==y&(cnt.lt$Age %in% 0:14)])/sum(cnt.lt$Lx[cnt.lt$Year==y&(cnt.lt$Age %in% 0:14)]),
                sum(cnt.lt$dx[cnt.lt$Year==y&(cnt.lt$Age %in% 15:64)])/sum(cnt.lt$Lx[cnt.lt$Year==y&(cnt.lt$Age %in% 15:64)]),
                sum(cnt.lt$dx[cnt.lt$Year==y&(cnt.lt$Age %in% 65:74)])/sum(cnt.lt$Lx[cnt.lt$Year==y&(cnt.lt$Age %in% 65:74)]),
                sum(cnt.lt$dx[cnt.lt$Year==y&(cnt.lt$Age %in% 75:84)])/sum(cnt.lt$Lx[cnt.lt$Year==y&(cnt.lt$Age %in% 75:84)]),
                sum(cnt.lt$dx[cnt.lt$Year==y&(cnt.lt$Age %in% 84:110)],na.rm=TRUE)/sum(cnt.lt$Lx[cnt.lt$Year==y&(cnt.lt$Age %in% 84:110)],na.rm=TRUE))
        hmd.ax <- nax(cnt.lt[which(cnt.lt$Year==min(y,max(cnt.lt$Year))),],age)
        e0.s <- get.lt(mx=mx.y,age=age,sex="M",hmd.ax)$ex[1]#life expectancy according to stmf data
        e0.h <- get.lt(mx=m.hmd,age=age,sex="M",hmd.ax)$ex[1]#life expectancy according to hmd data

        #ratios is the ratio between hmd rates and stmf. This is included as a quality check of stmf data
        ratios <- m.hmd/mx.y
    line <- c(cnt,y,e0.s,e0.h,ratios[1],ratios[2],ratios[3],ratios[4],ratios[5],ratios[6])
    mymatrix <- rbind(mymatrix,line)
    }
    #Compute 2020 life expectancy with two scenarios:
    # 1) mortality rates of future weeks are the same as those recorded on the same weeks in 2019
    # 2) mortality rates of future weeks are, the average of past weeks in 2020
    # Scenario 2 not really used in manuscript
    Av.L.Y <- max(as.numeric(names(table(stmf$Year[stmf$Country==cnt&stmf$Year<2020]))))
    Av.Weeks <- max(as.numeric(names(table(stmf$Week[stmf$Country==cnt&stmf$Year==2020]))))
    if(is.infinite(Av.Weeks)){Av.Weeks  <-  1}
    D.y.20 <- sum(stmf[which(stmf$CountryCode==cnt&stmf$Sex=="b"&stmf$Year==2020),10],na.rm=TRUE)
    D.nodate.20 <- input$Deaths[which(input$Sex=="b"&input$Year==2020&input$Week=="UNK"&input$Age=="TOT")]
    D.nodate.20 <- ifelse(length(D.nodate.20)==0,0,D.nodate.20)
    D.y.19 <- sum(stmf[which(stmf$CountryCode==cnt&stmf$Sex=="b"&stmf$Year==2019),10],na.rm=TRUE)
    D.nodate.19 <- input$Deaths[which(input$Sex=="b"&input$Year==2019&input$Week=="UNK"&input$Age=="TOT")]
    D.nodate.19 <- ifelse(length(D.nodate.19)==0,0,D.nodate.19)    
    if(cnt %in% c("FIN","NOR","USA","SWE")){Av.Weeks <- Av.Weeks-3}
    mx.2020 <- ((D.nodate.20+D.y.20)/D.y.20)*colSums(stmf[which(stmf$CountryCode==cnt&stmf$Sex=="b"&stmf$Year==2020&stmf$Week<=Av.Weeks),11:15])/Av.Weeks
    mx.2019 <- ((D.nodate.19+D.y.19)/D.y.19)*colSums(stmf[which(stmf$CountryCode==cnt&stmf$Sex=="b"&stmf$Year==Av.L.Y&stmf$Week>Av.Weeks),11:15])/max((52-Av.Weeks),1)
   mx.2020b <- ((D.nodate.20+D.y.20)/D.y.20)*colSums(stmf[which(stmf$CountryCode==cnt&stmf$Sex=="b"&stmf$Year==2020&stmf$Week%in%(10:(10+52-Av.Weeks))),11:15])/max((52-Av.Weeks),1)
    mh1 <- (Av.Weeks*mx.2020 + (52-Av.Weeks)*mx.2019)/52
    mh2 <- (Av.Weeks*mx.2020 + (52-Av.Weeks)*mx.2020b)/52
    e0.1 <- get.lt(mx=mh1,age=age,sex="M",hmd.ax)$ex[1]
    e0.2 <- get.lt(mx=mh2,age=age,sex="M",hmd.ax)$ex[1]
    line <- c(cnt,2020,e0.1,e0.2,NA,NA,NA,NA,NA,NA)
    mymatrix <- rbind(mymatrix,line)
}
mymatrix <- mymatrix[-1,]
mymatrix$e0.smtf[mymatrix$Country=="NZL_NP"&mymatrix$Year==2010] <- NA
