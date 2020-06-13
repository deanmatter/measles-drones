import numpy as np
np.set_printoptions(precision=3, suppress=True, linewidth=2000)
import random
from matplotlib import pyplot as plt
from matplotlib import patches
import copy
from collections import deque
import os
from json import loads, dumps

class vaccineDelivery:
    def __init__(self, d, eD):
        self.doses = d                  #the number of vaccine doses in this delivery
        self.expiryDate = eD            #the day t on which this delivery will expire

class POD:
    def __init__(self, name, cords, N, S, E, I, R, V, fT, mvpd, avt, mX, mY):
        self.name = name                #location at which POD is placed
        self.coordinates = cords        #tuple (x,y) of POD location in network
        self.N = N                      #total population 
        self.S = S                      #number of susceptible people
        self.E = E                      #number of exposed people
        self.last3E = deque()           #queue of the last 3 days' E values
        self.I = I                      #number of infectious people
        self.R = R                      #number of removed people
        self.vaccinated = V             #number of people vaccinated (by this intervention)
        self.deaths = 0                 #number of measles-induced deaths
        self.flightTime = fT            #number of minutes required to fly from DC to this POD and back
        self.vaccinesInStock = 0        #number of vaccine doses in stock at POD
        self.vaccsPerTeamDay = mvpd     #max number of vaccinations one team can do per day
        self.maxVaccinationsPerDay = 0  #this will be the total vaccs at this POD per day, for all teams
        self.vaccineDeliveries = deque()    #this is a queue data structure to store vaccine deliveries, FIFO.
        self.averageTurnout = avt       #the average population turnout for this POD for vaccination
        self.teamsAtPOD = 0
        self.mapX = mX
        self.mapY = mY
    
    def __str__(self):
        return self.name    
       
def readPODsFromFile(filename, maxVaccsTeamDay, turnout, targetedVaccination, vaccinationRate, populationMultiplier):
    '''Reads in location information from an input CSV file.'''
    PODs = []
    data = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), filename), 'r')
    lineCount = 0
    
    for line in data:
        if lineCount == 0:
            # skip the header row
            lineCount += 1
            continue
        lineCount += 1
        lineContents = line.strip().split(",")
        
        podName = lineContents[0]
        podTrnt = turnout
        podMaxVaccPD = maxVaccsTeamDay
        podmapX = int(lineContents[1])
        podmapY = int(lineContents[2])
        podXY = (podmapX,podmapY)                     #set x- and y- coords to unscaled values.

        podN = int(lineContents[3]) * populationMultiplier
        try: 
            podE = int(lineContents[4])
        except ValueError:
            # if there is a ValueError, set podE = 0
            podE = 0
        try:
            podI = int(lineContents[5])
        except ValueError:
            # if there is a ValueError, set podI = 0
            podI = 0

        #Split population according to vaccination rate
        podV, podS, podR = 0,0,0
        notInfected = podN - podE - podI
        if targetedVaccination:
            podV = int(notInfected * min(vaccinationRate,1))
            podS = notInfected - podV
        else:
            podR = int(notInfected * min(vaccinationRate,1))
            podS = notInfected - podR

        podflightTime = np.inf
        pod = POD(podName,podXY,podN,podS,podE,podI,podR,podV,podflightTime,podMaxVaccPD,podTrnt,podmapX,podmapY)
        PODs.append(pod)
    return PODs
    
def scaleCoordinatesAndDistances(unscaledDistances,PODs,maxDistance):
    ''' Returns: scaledDistances, and PODs with updated coordinates.'''
    #Calculate desired scaling factor and scale all distances accordingly.
    scaleFactor = maxDistance * 1.0 / np.amax(unscaledDistances)
    scaledDistances = unscaledDistances * scaleFactor
    
    #Scale the coordinates of each pod.
    for pod in PODs:
        scaledX = scaleFactor * pod.coordinates[0]
        scaledY = scaleFactor * pod.coordinates[1]
        pod.coordinates = (scaledX,scaledY)
      
    return scaledDistances, PODs

def findDC(PODs, podDistances):
    #maxInRow[i] stores the maximum distance from location i to any other location
    maxInRow = np.zeros(len(PODs))
    i = 0
    for pod in PODs:
        maxInRow[i] = max(podDistances[i])
        i += 1
    
    DCindex = -1   
    minimax = np.inf
    for i in range(0,len(maxInRow)):
        if maxInRow[i] < minimax:
            DCindex = i
            minimax = maxInRow[i]
    if minimax > 80:
        print("The best DC placement is at", PODs[DCindex], "but not all PODs are within 80km.")
    return DCindex
            
def calcPODdistanceMatrix(PODs):
    n = len(PODs)
    podDistances = np.zeros((n,n))
    i = 0
    for podi in PODs:
        j = 0
        ix, iy = podi.coordinates[0], podi.coordinates[1]
        for podj in PODs:
            jx, jy = podj.coordinates[0], podj.coordinates[1]
            podDistances[i][j] = np.sqrt((ix - jx)**2 + (iy - jy)**2)
            j += 1
        i += 1
    return podDistances
    
def calcPODflightTimes(PODs, podDistances, DCIndex, droneSpeed, flightLaunchTime):
    for i in range(0,len(PODs)):
        if i != DCIndex:
            PODs[i].flightTime = flightLaunchTime + (podDistances[DCIndex][i] * 2 / droneSpeed) * 60
        else:
            PODs[i].flightTime = 0

def calcMigrationFactor(PODs, i, j, podDistances, weights):
    if i == j:
        return 0
    alpha, beta, lamda, k = weights
    Tij = k * PODs[i].N**alpha * PODs[j].N**beta / podDistances[i][j]**lamda
    return Tij

def calcMigration(PODs, podDistances, migrationIntensity):
    '''Uses the gravity model of migration to calculate population proportions moving from POD i to j'''
    #weights = [0.25, 0.25, 2, migrationIntensity]       #-- Matadi
    weights = [0.5, 0.5, 2, migrationIntensity]        #--- Likasi
    migrations = np.zeros(np.shape(podDistances))
    for i in range(0,len(migrations)):
        for j in range(0,len(migrations[i])):
            migrations[i][j] = min(calcMigrationFactor(PODs, i, j, podDistances, weights), PODs[i].N)
        migrations[i] /= PODs[i].N
        migrations[i][i] = max(1 - sum(migrations[i]),0)
    #print(migrations)
    return migrations

def calcVaccinesNeeded(pod, numDays):
    unvaccinated = pod.S + pod.E + pod.R
    maxVaccsPossible = numDays * pod.maxVaccinationsPerDay
    turnout = pod.averageTurnout * numDays * pod.teamsAtPOD
    stock = pod.vaccinesInStock
    # Vaccines needed is number of likely vaccinations next3days, less the stock on hand.
    actualValue = max(0, min(unvaccinated, maxVaccsPossible, turnout) - stock)
    #multipleOf50 = np.ceil(actualValue/50) * 50     #since vaccines in boxes of 50
    #return int(multipleOf50)
    return int(actualValue)

def calcTurnout(pod):
    mean = pod.averageTurnout
    standDev = mean / 5
    turnoutVal = min(random.normalvariate(mean, standDev) * pod.teamsAtPOD, pod.N)
    return turnoutVal

def progressEpidemicByOneDay(PODs, params, MigrationProportions):
    n = len(PODs)
    Sarr = np.empty(n)
    Earr = np.empty(n)
    Iarr = np.empty(n)
    Rarr = np.empty(n)
    Varr = np.empty(n)
    deaths = 0
    for i in range(0,n):
        Sarr[i] = PODs[i].S
        Earr[i] = PODs[i].E
        Iarr[i] = PODs[i].I
        Rarr[i] = PODs[i].R
        Varr[i] = PODs[i].vaccinated
        deaths += PODs[i].deaths
    #print(Iarr)
          
    i = 0
    new_cases = np.zeros(len(PODs))
    for pod in PODs:
        #consider migrations
        migratePOD(pod, i, MigrationProportions, Sarr, Earr, Iarr, Rarr, Varr)
        #calculate SEIR progression for the day
        new_cases[i] = progressSinglePOD(pod, params)
        if len(pod.last3E) == 4:
            pod.last3E.popleft()
        i += 1 
    return PODs, new_cases
   
def progressSinglePOD(pod, params):
    beta, sigma, gamma, mu = params
    S = pod.S
    E = pod.E
    I = pod.I
    R = pod.R
    N = pod.N
    deaths = pod.deaths
    
    #calculation of changes in SEIR model
    newExposures = beta * S * I / N
    newInfectious = sigma * E
    newRecoveries = gamma * I
    newDeaths = mu * I
    pod.last3E.append(newExposures)
    
    #SEIR updates for the pod
    pod.S = S - newExposures 
    pod.E = E + newExposures - newInfectious 
    pod.I = I + newInfectious - newRecoveries - newDeaths
    pod.R = R + newRecoveries
    pod.deaths = deaths + newDeaths

    #Return new infectious as num of new cases today
    return newInfectious        
        
def migratePOD(pod, i, todaysMigration, Sarr, Earr, Iarr, Rarr, Varr):
    migrationIn = todaysMigration[:,i]
    pod.S = np.dot(Sarr, migrationIn)
    pod.E = np.dot(Earr, migrationIn)
    pod.I = np.dot(Iarr, migrationIn)
    pod.R = np.dot(Rarr, migrationIn)
    pod.vaccinated = np.dot(Varr, migrationIn)
    pod.N = pod.S + pod.E + pod.I + pod.R + pod.vaccinated + pod.deaths
      
def vaccinate(PODs):
    totVaccsGiven = 0
    vaccsGivenDC = 0
    
    vaccsGivenAtSpecificNode = 0    
    specificNodeName = ""         #just for checking a single node's vaccination progression
    
    for pod in PODs:
        vaccinesGiven = vaccinateOnePOD(pod)
        totVaccsGiven += vaccinesGiven
        if pod.vaccinesInStock == np.inf:
            vaccsGivenDC = vaccinesGiven
        #if vaccinesGiven > 0:
            #print(vaccinesGiven, "people vaccinated in", pod.name)
        if pod.name == specificNodeName:
            vaccsGivenAtSpecificNode = vaccinesGiven
    
    if specificNodeName != "":
        totVaccsGiven = vaccsGivenAtSpecificNode
            
    return totVaccsGiven, vaccsGivenDC

def vaccinateOnePOD(pod):    
    #print("\nVaccinating:-----------------------")
    if pod.maxVaccinationsPerDay == 0:
        return 0
    
    maxVaccinations = pod.maxVaccinationsPerDay
    
    podE72hrs = 0               #number of people in E that can be saved
    for i in range(0, len(pod.last3E)):
        podE72hrs += pod.last3E[i]
    
    turnout = calcTurnout(pod)
    if turnout <= 0:
        return 0
           
    #some vaccinations here are useless since they're given to people in R
    turnout = min(turnout, pod.S + pod.E + pod.R)
    if turnout == 0:
        return 0
    SturnoutProp = pod.S * 1.0 / (pod.S + pod.E + pod.R)
    EturnoutProp = pod.E * 1.0 / (pod.S + pod.E + pod.R)
    #print("Sprop",SturnoutProp,"Eprop",EturnoutProp)
    
    propE72hrs = 0              #proportion of E that is recently enough exposed to be vaccd
    if pod.E > 0:
        propE72hrs = (podE72hrs/pod.E)
    #print("E72/E prop",propE72hrs)
        
    vaccinesGiven = min(turnout, pod.vaccinesInStock, maxVaccinations)
    
    #print("vaccines Given",vaccinesGiven)

    numEffectiveSvaccinations = vaccinesGiven * SturnoutProp * vaccineEffectiveness
    #print("effective s",numEffectiveSvaccinations)
    numEffectiveEvaccinations = vaccinesGiven * EturnoutProp * propE72hrs *  prophylaxis72hrSuccessRate
    #print("effective E",numEffectiveEvaccinations)

    numRvaccinations = max(vaccinesGiven * (1-SturnoutProp-EturnoutProp),0)
    #print("number to R",numRvaccinations)
    pod.S -= numEffectiveSvaccinations
    pod.E -= numEffectiveEvaccinations
    pod.R -= numRvaccinations
    
    #pod.vaccinated += vaccinesGiven        -only add successfully vaccinated people to V class
    successfulVaccs = numEffectiveEvaccinations + numEffectiveSvaccinations + numRvaccinations
    pod.vaccinated += successfulVaccs
    #print("Successful vaccs",successfulVaccs)
    
    #Adjustment of E72hrs queue
    vaccinateE3(pod, numEffectiveEvaccinations)
    
    #Vaccine stock adjustment
    pod.vaccinesInStock -= vaccinesGiven    
    stillUnremoved = vaccinesGiven
    while stillUnremoved > 0 and pod.vaccinesInStock != np.inf and len(pod.vaccineDeliveries) > 0:
        if pod.vaccineDeliveries[0].doses > stillUnremoved:
            pod.vaccineDeliveries[0].doses -= stillUnremoved
            stillUnremoved = 0
        else: 
            stillUnremoved -= pod.vaccineDeliveries[0].doses
            pod.vaccineDeliveries.popleft()
            
    return vaccinesGiven

def vaccinateE3(pod, Evacc):
    sumE3 = sum(pod.last3E)
    if sumE3 == 0 or Evacc == 0:
        return
    l30vaccs = Evacc * (pod.last3E[0] / sumE3)  #vaccs given to people exposed yesterday
    l31vaccs = Evacc * (pod.last3E[1] / sumE3)  #vaccs given to people exposed 2 days ago
    pod.last3E[0] -= l30vaccs
    pod.last3E[1] -= l31vaccs
    pod.last3E[2] -= (Evacc - l30vaccs - l31vaccs)
    return

def expireVaccines(PODs, t):
    '''Assuming podDeliveries are all in FIFO order, removes vaccine cohorts that have expired.'''
    expiredVs = 0
    for pod in PODs:
        podDeliveries = pod.vaccineDeliveries
        if len(podDeliveries) == 0:
            continue
        
        while len(podDeliveries) > 0 and podDeliveries[0].expiryDate <= t:
            d = podDeliveries.popleft()
            pod.vaccinesInStock -= d.doses
            expiredVs += d.doses
            #print(d.doses, "vaccines have expired in", pod.name)
    return expiredVs

def deliverVaccinesByDroneSimple(vaccStrategy, t, PODs, workingMinutesPerDay, droneVC, numDrones, params, mdP):
    deliveries = 0
    vaccinesDelivered = 0
    times = np.zeros(numDrones)     #no launch staggering
    
    dronesFinished = np.full(numDrones, False)
    while all(dronesFinished) == False:
        for drone in range(0,numDrones):
            pod = choosePODtoFlyTo(vaccStrategy, PODs, times[drone], workingMinutesPerDay, mdP)
            if pod != -1:       # this ensures the pod needs vaccines, is not the DC, and can be delivered to in time
                times[drone] += pod.flightTime
                deliveryQty = droneVC
                pod.vaccinesInStock += deliveryQty
                deliveries += 1
                vaccinesDelivered += deliveryQty
                thisDelivery = vaccineDelivery(deliveryQty * 1, t + mdP)
                pod.vaccineDeliveries.append(thisDelivery)
                #if t in [100,101,102,103]:
                #    displayTime = str(int(times[drone]/60 + 7)) + ":" + "{:02d}".format(int(times[drone]%60))
                #    print(displayTime,"- Drone", drone+1, " delivered", deliveryQty,"vaccines to",
                #            pod.name, "which still needs", calcVaccinesNeeded(pod,mdP)
                #    )
            else:
                dronesFinished[drone] = True
    return deliveries, vaccinesDelivered

def deliverVaccinesByDroneEPE(t, PODs, workingMinutes, droneVC, numDrones, params, mdP):
    deliveries = 0
    vaccinesDelivered = 0
    
    n = len(PODs)
    ratios = np.zeros(n)
    for i in range(0,n):
        pod = PODs[i]
        if pod.vaccinesInStock != np.inf:
            ratios[i] = flightPreventedExposures(pod, droneVC, params, 1) / pod.flightTime

    time = 0
    droneAvail = np.zeros(numDrones)
    while time <= workingMinutes:
        #select POD with highest ratio of prevented exposures to flight time
        if np.max(ratios) > 0:
            podIndex = np.argmax(ratios)
            #ensure flight ends before the end of the day
            while PODs[podIndex].flightTime + time > workingMinutes:
                ratios[podIndex] = 0                # if this drone wouldn't return in time, set EPE/min=0
                podIndex = np.argmax(ratios)        # choose the second best EPE/min drone
                if max(ratios) == 0:                # if EPE/min=0 for all drones
                    time = workingMinutes + 1       # declare the day complete
                    break
        else:
            break
        pod = PODs[podIndex]
        
        deliveryQty = 0
        if calcVaccinesNeeded(pod,mdP) > 0:
            deliveryQty = droneVC

        if deliveryQty == 0:
            #the POD does not need any vaccines
            ratios[podIndex] = 0
            continue
    
        droneIndex = findAvailDrone(droneAvail, time, pod, workingMinutes)  #find available drone index -- returned, and will complete before workhours finish
        if droneIndex == -1:
            #no drone currently available
            time += 1       #check whether a drone will return the next minute
            continue
    
        #perform a drone delivery
        droneAvail[droneIndex] = time + pod.flightTime          #set the time that this drone will return
        deliveries += 1
        vaccinesDelivered += deliveryQty
        pod.vaccinesInStock += deliveryQty
        thisDelivery = vaccineDelivery(deliveryQty * 1, t + mdP)
        pod.vaccineDeliveries.append(thisDelivery)
        ratios[podIndex] = flightPreventedExposures(pod, droneVC, params, 1) / pod.flightTime
        #if t in [100,101,102,103]:
        #    displayTime = str(int((droneAvail[droneIndex])/60 + 7)) + ":" + "{:02d}".format((int(droneAvail[droneIndex]))%60)
        #    print(displayTime,"- Drone", droneIndex+1, " delivered", deliveryQty,
        #        "to", pod.name, "which now has",pod.vaccinesInStock,", which needs", calcVaccinesNeeded(pod,mdP), 
        #                    "vaccines. EPE updated to:", ratios[podIndex] * pod.flightTime          
        #        )
        #time += 1              #there is no need to increase (stagger) the time after a delivery anymore
    return deliveries, vaccinesDelivered
   
def findAvailDrone(droneAvail, time, pod, workMins):
    mI = np.argmin(droneAvail)  #select the drone that'll be (or has already) returning the soonest
    if droneAvail[mI] <= time and time + pod.flightTime <= workMins:    # if this drone has already returned and flying with it won't exceed working minutes
        return mI               #return the index of the drone
    else:
        return -1               #else return 'no drone available', essentially.
   
def deliverByDrone(vaccStrategy, t, PODs, workingMinutes, droneVC, numDrones, params, mdP):
    ''' Depending on vaccine delivery strategy, delivers vaccines by drone using helper methods.'''
    if vaccStrategy == 'EPE':     #clever exposure-prevention delivery schedule
        delivDetails = deliverVaccinesByDroneEPE(t, PODs, workingMinutes, droneVC, numDrones, params, mdP)
    elif vaccStrategy == 'uncapped':
        delivDetails = deliverUncappedVaccines(t, PODs, mdP)
    else:                   #simple delivery schedule
        delivDetails = deliverVaccinesByDroneSimple(vaccStrategy, t, PODs, workingMinutes, droneVC, numDrones, params, mdP)
    return delivDetails

def deliverUncappedVaccines(t, PODs, mdP):
    deliveries = 0
    vaccinesDelivered = 0
    deliveryQty = 999999                 # arbitrarily high. Can't be np.inf since that denotes the DC.
    for pod in PODs:
        deliveries += 1
        vaccinesDelivered += deliveryQty
        pod.vaccinesInStock += deliveryQty
        thisDelivery = vaccineDelivery(deliveryQty * 1, t + mdP)
        pod.vaccineDeliveries.append(thisDelivery)
    return deliveries, vaccinesDelivered

def choosePODtoFlyTo(vaccStrategy, PODs, time, workingMinutesPerDay, mdP):
    maxV = -1
    maxpod = -1
    for pod in PODs:
        #Different strategies for selecting the order of deliveries
        if pod.flightTime == 0:
            continue
        elif vaccStrategy == 'S':
            podVal = pod.S / pod.flightTime
        elif vaccStrategy == 'I':
            podVal = pod.I / pod.flightTime
        elif vaccStrategy == 'N':
            podVal = pod.N / pod.flightTime
        elif vaccStrategy == 'absS':
            podVal = pod.S
        elif vaccStrategy == 'absI':
            podVal = pod.I
        elif vaccStrategy == 'absN':
            podVal = pod.N
        else:
            print("Invalid delivery strategy selected. Choose absS, S, absI, I, absN, N, or EPE.")
            quit()
        
        #the drone flight must be able to return in time
        if time + pod.flightTime <= workingMinutesPerDay:
            #there's no need for more vaccines than people to vaccinate
            if calcVaccinesNeeded(pod, mdP) > 0:
                #no drone flights to DC
                if podVal > maxV and pod.flightTime != 0:
                    maxV = podVal
                    maxpod = pod
    
    if maxV == 0:
        return -1
    return maxpod
               
def flightPreventedExposures(pod, droneVC, params, numFlights):
    #create two variations of this POD, with differing vaccinesInStock values
    beforeFlightsPOD = copy.deepcopy(pod)
    afterFlightsPOD = copy.deepcopy(pod)
    afterFlightsPOD.vaccinesInStock += droneVC * numFlights
    thisDelivery = vaccineDelivery(droneVC * numFlights, np.inf)
    afterFlightsPOD.vaccineDeliveries.append(thisDelivery)
    
    #print("\nAfter delivery of", droneVC * numFlights)
    #print("Before",vars(beforeFlightsPOD))
    #print("After",vars(afterFlightsPOD))

    #vaccinate both PODs with the vaccines available at each
    vaccinateOnePOD(beforeFlightsPOD)
    vaccinateOnePOD(afterFlightsPOD)

    #print("\nAfter vaccination")
    #print("Before",vars(beforeFlightsPOD))
    #print("After",vars(afterFlightsPOD))
    
    #progress both PODs another day
    progressSinglePOD(beforeFlightsPOD, params)
    progressSinglePOD(afterFlightsPOD, params)

    #print("\nAfter progression")
    #print("Before",vars(beforeFlightsPOD))
    #print("After",vars(afterFlightsPOD))
    #print("Thus EPE",max(beforeFlightsPOD.E - afterFlightsPOD.E, 0))
    #print()
    #quit()
    
    #return the difference between exposure values - how many exposures did the flights prevent for the following day
    return max(beforeFlightsPOD.E - afterFlightsPOD.E, 0)
    
def teamPreventedExposures(pod, teams, params):
    #create two variations of this POD, with differing numbers of teams
    noTeamsPOD = copy.deepcopy(pod)
    afterTeamsPOD = copy.deepcopy(pod)
    
    noTeamsPOD.teamsAtPOD = 0
    noTeamsPOD.maxVaccinationsPerDay = 0
    
    afterTeamsPOD.teamsAtPOD = teams
    afterTeamsPOD.maxVaccinationsPerDay = afterTeamsPOD.vaccsPerTeamDay * afterTeamsPOD.teamsAtPOD
        
    #vaccinate both PODs with the vaccines available at each
    vaccinateOnePOD(noTeamsPOD)
    vaccinateOnePOD(afterTeamsPOD)
    
    #progress both PODs another day
    progressSinglePOD(noTeamsPOD, params)
    progressSinglePOD(afterTeamsPOD, params)
    
    #return the difference between exposure values - how many exposures did the flights prevent
    return max(noTeamsPOD.E - afterTeamsPOD.E, 0)
  
  
def assignTeams(PODs, numTeams, tStrategy):
    if tStrategy == 'spread':   
        # Splits teams equally among locations, randomly assigns remaining teams.
        remainingTeams = numTeams
        for pod in PODs:
            if remainingTeams > 0:
                pod.teamsAtPOD = min(roundUsingProb(numTeams/len(PODs)), remainingTeams)
                remainingTeams -= pod.teamsAtPOD
                pod.maxVaccinationsPerDay = pod.vaccsPerTeamDay * pod.teamsAtPOD
            else:
                break
            
        while remainingTeams > 0:
            randIndex = int(len(PODs) * random.random())
            PODs[randIndex].teamsAtPOD += 1
            remainingTeams -= 1
            PODs[randIndex].maxVaccinationsPerDay = PODs[randIndex].vaccsPerTeamDay * PODs[randIndex].teamsAtPOD
    elif tStrategy == 'S':    
        podSs = np.zeros(np.shape(PODs))
        for i in range(0,len(PODs)):
            PODs[i].teamsAtPOD = 0
            PODs[i].maxVaccinationsPerDay = 0
            podSs[i] = PODs[i].S
        teamsLeft = numTeams
        while teamsLeft > 0:
            bestSindex = np.argmax(podSs)
            podSproportion = podSs[bestSindex]*1.0 / sum(podSs)
            PODs[bestSindex].teamsAtPOD = roundUsingProb(teamsLeft * podSproportion)
            PODs[bestSindex].maxVaccinationsPerDay = PODs[bestSindex].vaccsPerTeamDay * PODs[bestSindex].teamsAtPOD
            teamsLeft -= PODs[bestSindex].teamsAtPOD
            #print(PODs[bestSindex].teamsAtPOD, "teams at", PODs[bestSindex].name)
            podSs[bestSindex] = 0            
    elif tStrategy == 'I':     
        podIs = np.zeros(np.shape(PODs))
        for i in range(0,len(PODs)):
            PODs[i].teamsAtPOD = 0
            PODs[i].maxVaccinationsPerDay = 0
            podIs[i] = PODs[i].I
        teamsLeft = numTeams
        while teamsLeft > 0:
            bestIindex = np.argmax(podIs)
            podIproportion = podIs[bestIindex]*1.0 / sum(podIs)
            PODs[bestIindex].teamsAtPOD = roundUsingProb(teamsLeft * podIproportion)
            PODs[bestIindex].maxVaccinationsPerDay = PODs[bestIindex].vaccsPerTeamDay * PODs[bestIindex].teamsAtPOD
            teamsLeft -= PODs[bestIindex].teamsAtPOD
            #print(PODs[bestIindex].teamsAtPOD, "teams at", PODs[bestIindex].name)
            podIs[bestIindex] = 0
    elif tStrategy == 'N':   
        podNs = np.zeros(np.shape(PODs))
        for i in range(0,len(PODs)):
            PODs[i].teamsAtPOD = 0
            PODs[i].maxVaccinationsPerDay = 0
            podNs[i] = PODs[i].N
        teamsLeft = numTeams
        while teamsLeft > 0:
            bestNindex = np.argmax(podNs)
            podNproportion = podNs[bestNindex]*1.0 / sum(podNs)
            PODs[bestNindex].teamsAtPOD = roundUsingProb(teamsLeft * podNproportion)
            PODs[bestNindex].maxVaccinationsPerDay = PODs[bestNindex].vaccsPerTeamDay * PODs[bestNindex].teamsAtPOD
            teamsLeft -= PODs[bestNindex].teamsAtPOD
            #print(PODs[bestNindex].teamsAtPOD, "teams at", PODs[bestNindex].name)
            podNs[bestNindex] = 0
    elif tStrategy == 'I/N':    
        podINs = np.zeros(np.shape(PODs))
        for i in range(0,len(PODs)):
            PODs[i].teamsAtPOD = 0
            PODs[i].maxVaccinationsPerDay = 0
            podINs[i] = PODs[i].I * 1.0 / PODs[i].N
        teamsLeft = numTeams
        while teamsLeft > 0:
            bestINindex = np.argmax(podINs)
            podINproportion = podINs[bestINindex]*1.0 / sum(podINs)
            PODs[bestINindex].teamsAtPOD = roundUsingProb(teamsLeft * podINproportion)
            PODs[bestINindex].maxVaccinationsPerDay = PODs[bestINindex].vaccsPerTeamDay * PODs[bestINindex].teamsAtPOD
            teamsLeft -= PODs[bestINindex].teamsAtPOD
            #print(PODs[bestINindex].teamsAtPOD, "teams at", PODs[bestINindex].name)
            podINs[bestINindex] = 0
            
    elif tStrategy == 'EPE':   
        podEPEs = np.zeros((len(PODs),numTeams+1))
        podEPEperTeam = np.zeros((len(PODs),numTeams+1))
        podTeamsAssigned = np.zeros(len(PODs),dtype=int)    #current teams assigned to each pod
        add1teamEPEs = np.zeros(len(PODs))    #array of EPE vals for adding one more team at each pod
        for nT in range(1,numTeams+1):        #build each column before the row
            for podI in range(0,len(PODs)):
                if nT == 1 or podEPEs[podI][1] != 0:        #don't calculate this excessively
                    podEPEs[podI][nT] = teamPreventedExposures(PODs[podI], nT, params)
                
                if podEPEs[podI][1] != 0:      #if adding one team doesn't prevent an exposure, don't add more.
                    podEPEperTeam[podI][nT] = max(podEPEs[podI][nT] - podEPEs[podI][nT-1], 0)
                
                teamsAssignedToPODalready = podTeamsAssigned[podI]
                add1teamEPEs[podI] = podEPEperTeam[podI][teamsAssignedToPODalready+1]
            newTeamPODindex = np.argmax(add1teamEPEs)
            podTeamsAssigned[newTeamPODindex] += 1
            
            if sum(podTeamsAssigned) == numTeams:
                break    
        
        for i in range(0,len(PODs)):
            PODs[i].teamsAtPOD = podTeamsAssigned[i]
            PODs[i].maxVaccinationsPerDay = PODs[i].vaccsPerTeamDay * podTeamsAssigned[i]
            
    else:
        print("Invalid selection of team strategy")
        quit()
  
def updatePlotHistory(PODs, plots):
    i = 0
    for pod in PODs:
        plots[i][0].append(pod.S)
        plots[i][1].append(pod.E)
        plots[i][2].append(pod.I)
        plots[i][3].append(pod.R)
        plots[i][4].append(pod.deaths)
        plots[i][5].append(pod.vaccinated)
        i += 1
 
def roundUsingProb(val):
    '''Rounds a float up or down randomly and probabilistically, depending on
    how close the value is to its floor or ceiling.'''
    if val % 1 > random.random():
        val = int(val) + 1
    else:
        val = int(val)
    return val
   
def plotPOD(daysOfIntervention, plots, podIndex, podName):
    plt.plot(np.arange(daysOfIntervention), plots[podIndex][0], label='Susceptible')
    plt.plot(np.arange(daysOfIntervention), plots[podIndex][1], label='Exposed')
    plt.plot(np.arange(daysOfIntervention), plots[podIndex][2], label='Infectious')
    #plt.plot(np.arange(daysOfIntervention), plots[podIndex][3], label='Recovered')
    plt.plot(np.arange(daysOfIntervention), plots[podIndex][4], label='Deaths')
    plt.plot(np.arange(daysOfIntervention), (plots[podIndex][5]-min(plots[podIndex][5])), label='Vaccinations')
    plt.title(label=podName)
    plt.ylabel("Number of People")
    plt.xlabel("Time Elapsed (Days)")
    plt.legend()
    plt.show()
    
def plotPODSum(daysOfIntervention, plots, podIndexes, PODs):
    Ss = np.array(plots[podIndexes[0]][0])
    Es = np.array(plots[podIndexes[0]][1])
    Is = np.array(plots[podIndexes[0]][2])
    Rs = np.array(plots[podIndexes[0]][3])
    deads = np.array(plots[podIndexes[0]][4])
    vacs = np.array(plots[podIndexes[0]][5])
    for i in range(1, len(podIndexes)):
        Ss = Ss + plots[podIndexes[i]][0]
        Es = Es + plots[podIndexes[i]][1]
        Is = Is + plots[podIndexes[i]][2]
        Rs = Rs + plots[podIndexes[i]][3]
        deads = deads + plots[podIndexes[i]][4]
        vacs = vacs + plots[podIndexes[i]][5]
        
    plt.plot(np.arange(daysOfIntervention), Ss, label='Susceptible',color='gold')
    plt.plot(np.arange(daysOfIntervention), Es, label='Exposed',color='darkorange')
    plt.plot(np.arange(daysOfIntervention), Is, label='Infectious',color='crimson')
    plt.plot(np.arange(daysOfIntervention), deads, label='Deaths',color='k')
    plt.plot(np.arange(daysOfIntervention), Rs, label='Recovered',color='forestgreen')

    # Only plot the vaccinations curve if some where performed
    if max(vacs) > 0:
        plt.plot(np.arange(daysOfIntervention), (vacs-min(vacs)), label='Vaccinations',color='cornflowerblue')

    label = ""
    for i in podIndexes:
        label = label + PODs[i].name + ","
    label = label[:len(label)-1]
    if len(podIndexes) == len(PODs):
        label = "Measles epidemic progression across the whole network"
    
    plt.title(label)
    #plt.title("Progression of measles epidemic SEIR model")
    plt.ylabel("Number of people")
    plt.xlabel("Day of epidemic")
    plt.legend(loc='right')
    
    import matplotlib as mpl
    mpl.rcParams['figure.dpi'] = 150
    plt.savefig(f"{plot_filename}.pdf",bbox_inches='tight')
    #plt.show()
    plt.close()

    
def plotMap(PODs, t, waitingForIntervention, IstartTime, intOver, 
            totExpired, totVaccs, totVaccsGiven, totCost):
    image = plt.imread("map.png")
    fig, ax = plt.subplots(figsize=(12,10))
    ax.imshow(image)
    ax.axis('off')
    for pod in PODs:
        if pod.I == 0:
            ax.text((pod.mapX), (pod.mapY + 40), pod.name, fontsize = 12, fontweight='normal', horizontalalignment='center')
        elif pod.vaccinesInStock > 0 and intOver == False:
            ax.text((pod.mapX), (pod.mapY + 40), pod.name, fontsize = 12, fontweight='normal', horizontalalignment='center', color='seagreen')
        else:
            ax.text((pod.mapX), (pod.mapY + 40), pod.name, fontsize = 12, fontweight='normal', horizontalalignment='center', color='darkred')
        
        vaccRate = int(pod.vaccinated / (pod.N-pod.deaths) * 100)
        ax.text((pod.mapX), (pod.mapY + 80), str(vaccRate) + "%", fontsize = 10, fontweight='light', horizontalalignment='center', color='seagreen')
        
        if intOver == False:
            for i in range(0,pod.teamsAtPOD):
                #plots the number of teams at the POD
                ax.plot(pod.mapX-45 - np.floor(i/5)*20, (pod.mapY - (i%5)*20), "o", markersize = 4, color='k')
            
        #plots the bar of the number of deaths at the POD
        if pod.I > 0:
            infectRect = patches.Rectangle((pod.mapX+5,pod.mapY), 20, pod.I/30, angle=180, color='goldenrod')
            ax.add_patch(infectRect)
            
        #plots the bar of the number of deaths at the POD
        if pod.deaths > 0:
            deathRect = patches.Rectangle((pod.mapX+45,pod.mapY), 20, pod.deaths/30, angle=180, color='darkred')
            ax.add_patch(deathRect)
            
    totDeaths = totI = 0
    for pod in PODs:
        totDeaths += pod.deaths
        totI += pod.I
    
    #show numbers for vaccinations, infectious and deaths, and day
    if totVaccs > 0:
        ax.text(ax.get_xlim()[1]-20, ax.get_ylim()[1]-90,"Vaccinations: " + str(totVaccsGiven), horizontalalignment='right', fontsize=17, fontweight='light', color='seagreen')
        ax.text(ax.get_xlim()[1]-20, ax.get_ylim()[1]-40,"(" + str(int(1000*totExpired/totVaccs)/10) + "% wastage)", horizontalalignment='right', fontsize=15, fontweight='light', color='seagreen')
    else:
        ax.text(ax.get_xlim()[1]-20, ax.get_ylim()[1]-40,"Vaccinations: " + str(totVaccsGiven), horizontalalignment='right', fontsize=17, fontweight='light', color='seagreen')
    ax.text(ax.get_xlim()[1]-20, ax.get_ylim()[1]+30,"Infectious: " + str(totI),horizontalalignment='right', fontsize=17, fontweight='light', color='goldenrod')
    ax.text(ax.get_xlim()[1]-20, ax.get_ylim()[1]+100,"Deaths: " + str(totDeaths),horizontalalignment='right', fontsize=17, fontweight='light', color='darkred')
    ax.text(ax.get_xlim()[1]-20, ax.get_ylim()[0],"Day: " + str(t),horizontalalignment='right', fontsize=18, fontweight='light', color='k')
    
    #show cost details and simulation status
    if waitingForIntervention == True:
        #waiting for intervention on day IstartTime
        ax.text(ax.get_xlim()[0], ax.get_ylim()[1]-90,"Days until intervention: " + str(IstartTime-t),horizontalalignment='left', fontsize=17, fontweight='light')
    elif intOver == True:
        ax.text(ax.get_xlim()[0], ax.get_ylim()[1]-90,"Intervention Complete", horizontalalignment='left', fontsize=17, fontweight='light')
        
    ax.text(ax.get_xlim()[0], ax.get_ylim()[1]-20,"Total Cost: $" + str(int(totCost)), horizontalalignment='left', fontsize=17, fontweight='light')

    print("saving figure", t)
    plt.savefig("animation/"+str(t)+".png")
    #plt.show()
    plt.close()
    

def simulate(filename='Likasi.csv'):
    #===============================================================================
    # Initial Calculations
    #===============================================================================
    PODs = readPODsFromFile(filename, maxVaccsTeamDay, turnout, targetedVaccination, vaccinationRate, populationMultiplier)
    unscaledDistances = calcPODdistanceMatrix(PODs)
    podDistances, PODs = scaleCoordinatesAndDistances(unscaledDistances, PODs, maxDistance)
    MigrationProportions = calcMigration(PODs, podDistances, migrationIntensity)
    DCindex = findDC(PODs, podDistances)
    #print("The DC is placed in", PODs[DCindex])
    PODs[DCindex].vaccinesInStock = np.inf
    calcPODflightTimes(PODs, podDistances, DCindex, droneSpeed, flightLaunchTime)
    waitingForIntervention = False
    interventionOver = False
    interventionStartTime = np.inf
    plots = []
    for i in range(0,len(PODs)):
        plots.append([[],[],[],[],[],[]])
    totDroneDelvs = 0                            #total number of drone flights for simulation
    totVaccs = 0                            #total number of vaccines delivered, incl DC vaccines
    totExpired = 0                          #total number of vaccines expired    
    totVaccsGiven = 0                       #total number of vaccines actually administered, incl DC.
    vaccsPerDay = np.zeros(simulationRuntime)
    cumulative_cases_pods = np.zeros(len(PODs))   # number of infectious cases so far, per pod

    #===============================================================================
    # Simulation - Main Daily Loop
    #===============================================================================    
    for t in range(0, simulationRuntime):
        #print("\n",t,":")
            
        #Progress the SEIR models and migration by one day
        PODs, new_cases = progressEpidemicByOneDay(PODs, params, MigrationProportions)
        cumulative_cases_pods += new_cases

        #if the intervention is running today
        if t >= interventionStartTime and t < interventionStartTime + interventionLength and random.random() < workDaysPerWeek/7:
            waitingForIntervention = False
            
            #expire old vaccines
            totExpired += expireVaccines(PODs, t)
            
            #assign team locations
            assignTeams(PODs, numTeams, teamStrategy)
            
            #deliver vaccines
            if deliveryType == "drone":
                delivDetails = deliverByDrone(vaccStrategy, t, PODs, workingMinutesPerDay, 
                                            droneVaccineCapacity, numberOfDrones, params, 
                                            monoDaysPotency)
                totDroneDelvs += delivDetails[0]
                totVaccs += delivDetails[1]
            elif deliveryType != "none":            #no vaccinations happen if deliveryType is none.
                print("Invalid delivery type selected.")
                quit()        
            
            #process vaccinations for the day
            if deliveryType != "none":
                vGiven, vGivenDC = vaccinate(PODs)
                vaccsPerDay[t] = vGiven
                totVaccsGiven += vGiven
                totVaccs += vGivenDC        # Add the vaccs given at DC to the total cumul. vaccines delivered
        elif t >= interventionStartTime + interventionLength:
            #intervention is over
            interventionOver = True
            #print("Intervention is over")
        elif not waitingForIntervention and t <= interventionStartTime:
            for pod in PODs:
                if pod.I/pod.N > interventionCaseRatio:
                    #if a confirmed measles case found in a town, then decide to intervene.
                    interventionStartTime = t + interventionLeadTime
                    waitingForIntervention = True
                    #print("The epidemic has been detected in", pod.name, "! Intervention will begin on day", interventionStartTime)            
                    break
        
        deliveryCost = totDroneDelvs * costPerFlight
        vaccineCost = totVaccs * costPerDoseMono
            
        updatePlotHistory(PODs, plots)
        
        #plot this iteration's png
        #plotMap(PODs, t, waitingForIntervention, interventionStartTime, interventionOver, 
        #       totExpired, totVaccs, totVaccsGiven, deliveryCost + vaccineCost)
    
    # ----------------- Result Reporting
    deaths = 0
    for pod in PODs:
        deaths += pod.deaths
    #print("\nThe total number of deaths is:", deaths)
    
    # The total cases reported is the number of infectious people simulated
    total_cases = sum(cumulative_cases_pods)

    #print("Total cost of monodose vaccines delivered: $", vaccineCost, ",for", totVaccs, "doses.")
    if totVaccs > 0:
        expiryRatio = totExpired / totVaccs * 100
        #print("Percentage of vaccines expired without use:", round(expiryRatio,2), "%")

    #print("Total number of vaccines delivered, plus those used at the DC:",totVaccs)
    #print("Total number of vaccines actually administered to patients:", totVaccsGiven)
    
    #plotPODSum(simulationRuntime, plots, np.arange(0,len(PODs)), PODs)
    #plotPOD(simulationRuntime, plots, 7, "Likasi")
    #plotMap(PODs,  t, waitingForIntervention, interventionStartTime, interventionOver,
    #      totExpired, totVaccs, totVaccsGiven, vaccineCost + deliveryCost)
    
    if vaccStrategy == 'uncapped':
        # In this case, the total vaccines delivered and used is just the total vaccs given (which includes DC)
        totVaccs = totVaccsGiven

    return total_cases, deaths, totVaccs, deliveryCost + vaccineCost  #totExpired, vaccsPerDay #,Vactots,Stots


def simulateRepeatedly(filename, repetitions=50):
    ''' Performs repeated simulations, returning the averages of the result metrics. '''
    ca, d, v, co = (np.zeros(repetitions) for i in range(4))    # creates one np zero array each
    for i in range(repetitions):
        ca[i],d[i],v[i],co[i] = simulate(filename)
    
    # Check that standard deviation is acceptable
    for arr in [ca,d,v,co]:
        mean = np.average(arr)
        stdev = np.std(arr,ddof=1)          # 1 degree of freedom for sample stdev
        Zval = 1.96                         # using alpha = 0.05
        error = 0.05                        # denoted by epsilon in formula

        if mean == 0:
            continue

        sims_required = ((Zval*stdev)/(error*mean))**2
        if sims_required > repetitions:
            print(f"ERROR: Insufficient simulations:{sims_required} needed, {repetitions} done.")
            return(-1,-1,-1,-1)
    return np.average(ca),np.average(d),np.average(v),np.average(co)
    

#Parameters    ========================================================================
simulationRuntime = 280             #days to run the simulation for
#Measles SEIR parameters
exposedDays = 10                    #number of days a patient is exposed for without symptoms
infectiousDays = 8                  #number of days a patient is infectious for
deathRate = 0.0329 * 1/infectiousDays #daily death rate. From wolfson2009estimates - 3.29 mean, 0.05-6% WHO estimate for low income countries
R0 = 15                             #basic reproductive number of the epidemic
params = [R0 / infectiousDays, 1/exposedDays, 1/infectiousDays, deathRate]
migrationIntensity = 1              #factor by which migration is multiplied. 2 means more migration.
#vaccine parameters
vaccineEffectiveness = 0.95         #probability the vaccine works (for non-exposed)
prophylaxis72hrSuccessRate = 0.83   #probability the vaccine works (for exposed, within 72hrs)
monoDaysPotency = 3                 #number of days for which the vaccine lasts outside of cold-chain
#intervention parameters
interventionLeadTime = 15           #number of days before vaccination starts
interventionCaseRatio = 0.009       #ratio of I/S in a town before detection
interventionLength = np.inf         #number of days the intervention lasts for
workingMinutesPerDay = 660          #11 working hours per day: 7am to 6pm
workDaysPerWeek = 7                 #number of working days per week for MSF teams
numTeams = 15                       #number of vaccination teams in the field
maxVaccsTeamDay = 2000              #teams can vaccinate up to ~2000 per day - poncin2018implementation
turnout = 999999                    #turnout was 900, now inf to effectively remove its impact.
#delivery details
flightLaunchTime = 10               #minutes per flight, to set up takeoff
droneSpeed = 100                    #100 kilometres per hour     
numberOfDrones = 5                  #number of drones
droneVaccineCapacity = 60           #number of vaccine doses per drone
#costs
costPerDoseMono = 2.85              #the cost per dose of monodose measles vaccine
costPerFlight = 17                  #$17 per drone flight
#strategies
vaccStrategy = 'N'                  #I, S, N, EPE, uncapped, absI, absS, absN
teamStrategy = 'N'                  #I, S, N, EPE, I/N, spread
deliveryType = 'drone'               #"none", "drone"
targetedVaccination = False         #True: already-vaccd people go to V. False: they go to R category.
#input dataset
maxDistance = 40                    #The distance in km that the max inter-location distance is scaled to
vaccinationRate = 0.66              #66% vaccination rate in network
populationMultiplier = 1            #Used for sensitivity analysis only. Multiplier for population size

# Sensitivity Analysis
epidemic_parameters = [
    "exposedDays",
    "infectiousDays",
    "R0",
    "deathRate"
]

vaccination_parameters = [
    "vaccineEffectiveness",
    "prophylaxis72hrSuccessRate",
    "monoDaysPotency",
    "vaccinationRate"
]

intervention_parameters = [
    "numTeams",
    "workingMinutesPerDay",
    "maxVaccsTeamDay",
    "interventionCaseRatio",
    "interventionLeadTime"
]

network_parameters = [
    "migrationIntensity",
    "maxDistance",
    "populationMultiplier"
]

drone_parameters = [
    "droneVaccineCapacity",
    "numberOfDrones",
    "droneSpeed",
    "flightLaunchTime"
]

network_type = 'monocentric'
maxDistance = 40
output_folder = f"results/sensitivity"

with open(f"{output_folder}/{network_type}.csv","a+") as f:
    f.write("Parameter,-50%,,,,-30%,,,,-20%,,,,-10%,,,,+10%,,,,+20%,,,,+30%,,,,+50%,,,,")
    for param_set in [epidemic_parameters,vaccination_parameters,intervention_parameters,network_parameters,drone_parameters]:
        f.write("\n")
        for p_name in param_set:
            # Print parameter name to file
            f.write(f"{p_name},")

            # Set p_val = original variable value
            p_val = vars()[p_name]
            print(f"Testing {p_name}, with initial value {p_val}")

            for i in [0.5*p_val,0.7*p_val, 0.8*p_val, 0.9*p_val, 1.1*p_val, 1.2*p_val, 1.3*p_val,1.5*p_val]:
                # Set the variable's value to i for each iteration, potentially with truncation to integer
                if p_name in ['numTeams','interventionLeadTime','numberOfDrones']:
                    i = int(i)
                vars()[p_name] = i
                print(f"Set {p_name}={i}")

                # Update the parameters to be used
                params = [R0 / infectiousDays, 1/exposedDays, 1/infectiousDays, deathRate]
                
                #Simulate this epidemic and write results to file
                c,d,v,cost = simulateRepeatedly(f"Generic_network_{network_type}.csv", 1)
                f.write(f"{c},{d},{v},{cost},")
            vars()[p_name] = p_val
            print(f"Reset {p_name}={p_val}")
            f.write("\n")
        
#TODO: (Optional) Randomly generated networks
#TODO: (Optional) Lockdown scenario of less migration between nodes. Reduced Ro?
