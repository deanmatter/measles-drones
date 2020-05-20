import numpy as np
np.set_printoptions(precision=1, suppress=True, linewidth=2000)
import random
from matplotlib import pyplot as plt
from matplotlib import patches
import copy
from collections import deque

class vaccineDelivery:
    def __init__(self, d, eD):
        self.doses = d                  #the number of vaccine doses in this delivery
        self.expiryDate = eD            #the day t on which this delivery will expire

class POD:
    def __init__(self, name, urban, cords, N, S, E, I, R, V, fT, mvpd, avt, mX, mY):
        self.name = name                #location at which POD is placed
        self.isUrban = urban           #true/false value of if the DC can be placed at this POD
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
       
def readPODsFromFile(filename, mvpFreeDay, mvpFixedDay, ruralT, urbanT):
    '''Reads in location information from an input CSV file.'''
    PODs = []
    roadDistances = []
    data = open(filename,'r')
    lineCount = -2
    mode = 0
    
    for line in data:
        lineCount += 1
        lineContents = line.strip().split(",")
        
        if lineCount == -1:
            mode = 1
            continue
        if lineContents[0][:1] == "#":
            mode += 1
            continue
        
        #Mode 1 -- reading in POD details
        if mode == 1:
            podName = lineContents[0]
            podUrban = False
            podTrnt = ruralT
            podMaxVaccPD = mvpFreeDay
            if int(lineContents[1]) == 1:
                podUrban = True
                podTrnt = urbanT
                podMaxVaccPD = mvpFixedDay
            podXY = (float(lineContents[2]),float(lineContents[3]))
            podN = int(lineContents[4])
            podS = int(lineContents[5])
            podE = int(lineContents[6])
            podI = int(lineContents[7])
            podR = int(lineContents[8])
            podV = int(lineContents[9])
            podmapX = int(lineContents[10])
            podmapY = int(lineContents[11])
            podfT = np.inf
            pod = POD(podName,podUrban,podXY,podN,podS,podE,podI,podR,podV,podfT,podMaxVaccPD,podTrnt,podmapX,podmapY)
            PODs.append(pod)
                
        #Mode 2 -- reading in road times details
        if mode == 2 and lineContents[0] != '':
            thisRowsRoadDs = []
            for i in range(1,len(PODs)+1):
                if(lineContents[i] == ''):
                    thisRowsRoadDs.append(np.inf)
                else:
                    thisRowsRoadDs.append(float(lineContents[i]))
            roadDistances.append(thisRowsRoadDs)
            
    roadDistances = np.array(roadDistances)
    return PODs, roadDistances
    
    
def findDC(PODs, podDistances):
    #maxInRow[i] stores the maximum distance from location i to any other location
    maxInRow = np.zeros(len(PODs))
    i = 0
    for pod in PODs:
        if pod.isUrban == True:
            maxInRow[i] = max(podDistances[i])
        else:
            maxInRow[i] = np.inf
        i += 1
    
    DCindex = -1   
    minimax = np.inf
    for i in range(0,len(maxInRow)):
        if maxInRow[i] < minimax:
            DCindex = i
            minimax = maxInRow[i]
    if minimax > 80:
        print("The best DC placement is at", PODs[DCindex], "but not all PODs are within 80km!")
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

def calcRoadTimes(PODs, roadDistances, podDistances, roadCloseFactor):
    '''Returns two nxn matrices, with i-j entries =1 if the road i-j is open, and =0 if the road is closed. Also, road times'''
    roadTimes = np.zeros(np.shape(roadDistances))
    openRoads = np.zeros(np.shape(roadDistances))
    
    numClosedRoads = 0
    for r in range(0,len(PODs)):
        for c in range(0, len(PODs)):
            if roadDistances[r][c] == np.inf:
                #there is no road between r and c - 10kph travel speed by foot
                openRoads[r][c] = 0
            elif PODs[r].isUrban == True and PODs[c].isUrban == True:
                #travel speed is 60kph - between urban and urban
                roadTimes[r][c] = int(roadDistances[r][c] / 60 * 60)
                openRoads[r][c] = min(roundUsingProb(1 / roadCloseFactor),1)
            elif PODs[r].isUrban == False and PODs[c].isUrban == False:
                #travel speed is 20 kph - between rural and rural
                roadTimes[r][c] = int(roadDistances[r][c] / 20 * 60)
                openRoads[r][c] = min(roundUsingProb(0.6 / roadCloseFactor),1)
            else:
                #travel speed is 40kph - between rural and urban
                roadTimes[r][c] = int(roadDistances[r][c] / 40 * 60)
                openRoads[r][c] = min(roundUsingProb(0.8 / roadCloseFactor),1)
            
            if openRoads[r][c] == 0:
                roadTimes[r][c] = int(podDistances[r][c] / 10 * 60)
    
    for r in range(0,len(PODs)):
        for c in range(0, len(PODs)):
            if openRoads[r][c] == 0 and roadDistances[r][c] != np.inf:
                numClosedRoads += 1
            
    return numClosedRoads, openRoads, roadTimes

def calcMaxVaccinesNeeded50(pod, numDays):
    unvaccinated = pod.S + pod.E + pod.R
    maxVaccsPossible = numDays * pod.maxVaccinationsPerDay
    turnout = pod.averageTurnout * numDays * pod.teamsAtPOD
    stock = pod.vaccinesInStock
    actualValue = max(0, min(unvaccinated, maxVaccsPossible, turnout) - stock)
    #multipleOf50 = np.ceil(actualValue/50) * 50     #since vaccines in boxes of 50
    multipleOf50 = actualValue
    return int(multipleOf50)

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
    #print(sum(Iarr))
          
    i = 0
    for pod in PODs:
        #consider migrations
        migratePOD(pod, i, MigrationProportions, Sarr, Earr, Iarr, Rarr, Varr)
        #calculate SEIR progression for the day
        progressSinglePOD(pod, params)
        if len(pod.last3E) == 4:
            pod.last3E.popleft()
        i += 1
    return PODs
   
def progressSinglePOD(pod, params):
    beta, sigma, gamma, mu = params
    S = pod.S
    E = pod.E
    I = pod.I
    R = pod.R
    N = pod.N
    deaths = pod.deaths
    
    #calculation of changes in SEIR model
    newExposures = roundUsingProb(beta * S * I / N)
    newInfectious = roundUsingProb(sigma * E)
    newRecoveries = roundUsingProb(gamma * I)
    newDeaths = roundUsingProb(mu * I)
    pod.last3E.append(newExposures)
    
    #SEIR updates for the pod
    pod.S = S - newExposures 
    pod.E = E + newExposures - newInfectious 
    pod.I = I + newInfectious - newRecoveries - newDeaths
    pod.R = R + newRecoveries
    pod.deaths = deaths + newDeaths
        
def migratePOD(pod, i, todaysMigration, Sarr, Earr, Iarr, Rarr, Varr):
    migrationIn = todaysMigration[:,i]
    pod.S = roundUsingProb(np.dot(Sarr, migrationIn))
    pod.E = roundUsingProb(np.dot(Earr, migrationIn))
    pod.I = roundUsingProb(np.dot(Iarr, migrationIn))
    pod.R = roundUsingProb(np.dot(Rarr, migrationIn))
    pod.vaccinated = roundUsingProb(np.dot(Varr, migrationIn))
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
    
    propE72hrs = 0              #proportion of E that is recently enough exposed to be vaccd
    if pod.E > 0:
        propE72hrs = (podE72hrs/pod.E)
        
    vaccinesGiven = int(min(turnout, pod.vaccinesInStock, maxVaccinations))      
    
    numEffectiveSvaccinations = int(vaccinesGiven * SturnoutProp * vaccineEffectiveness)
    numEffectiveEvaccinations = int(vaccinesGiven * EturnoutProp * propE72hrs *  prophylaxis72hrSuccessRate)
    #numRvaccinations = vaccinesGiven - numEffectiveEvaccinations - numEffectiveSvaccinations
    numRvaccinations = max(int(vaccinesGiven * (1-SturnoutProp-EturnoutProp)),0)
    pod.S -= numEffectiveSvaccinations
    pod.E -= numEffectiveEvaccinations
    pod.R -= numRvaccinations
    
    #pod.vaccinated += vaccinesGiven        -only add successfully vaccinated people to V class
    successfulVaccs = numEffectiveEvaccinations + numEffectiveSvaccinations + numRvaccinations
    pod.vaccinated += successfulVaccs
    
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
    l30vaccs = roundUsingProb(Evacc * (pod.last3E[0] / sumE3))  #vaccs given to people exposed yesterday
    l31vaccs = roundUsingProb(Evacc * (pod.last3E[1] / sumE3))  #vaccs given to people exposed 2 days ago
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

def deliverVaccinesByDroneSimple(strategy, t, PODs, workingMinutesPerDay, droneVC, numDrones, params, mdP):
    deliveries = 0
    vaccinesDelivered = 0
    times = np.zeros(numDrones)
    #launch drones in staggered 5-minute intervals
    for i in range(1,len(times)):
        times[i] = times[i-1] + 5
    
    dronesFinished = np.full(numDrones, False)
    while all(dronesFinished) == False:
        for drone in range(0,numDrones):
            pod = choosePODtoFlyTo(strategy, PODs, times[drone], workingMinutesPerDay, mdP)
            if pod != -1:
                times[drone] += pod.flightTime
                deliveryQty = min(60*np.ceil(calcMaxVaccinesNeeded50(pod,mdP)/60), droneVC)
                if deliveryQty == 0:
                    continue
                pod.vaccinesInStock += deliveryQty
                deliveries += 1
                vaccinesDelivered += deliveryQty
                thisDelivery = vaccineDelivery(deliveryQty * 1, t + mdP)
                pod.vaccineDeliveries.append(thisDelivery)
                displayTime = str(int(times[drone]/60 + 7)) + ":" + "{:02d}".format(int(times[drone]%60))
                #print(displayTime,"- Drone", drone+1, " delivered", deliveryQty,"vaccines to", pod.name)
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
        if np.max(ratios) > 0.01:
            podIndex = np.argmax(ratios)
            #ensure flight ends before the end of the day
            while PODs[podIndex].flightTime + time > workingMinutes:
                ratios[podIndex] = 0
                podIndex = np.argmax(ratios)
                if max(ratios) == 0:
                    time = workingMinutes + 1
                    break
        else:
            break
        pod = PODs[podIndex]
        
        deliveryQty = min(np.ceil(calcMaxVaccinesNeeded50(pod,mdP)/60)*60, droneVC)
        if deliveryQty == 0:
            #the POD does not need any vaccines
            ratios[podIndex] = 0
            continue
    
        droneIndex = findAvailDrone(droneAvail, time, pod, workingMinutes)
        if droneIndex == -1:
            #no drone currently available
            time += 5
            continue
    
        #perform a drone delivery
        droneAvail[droneIndex] = time + pod.flightTime
        deliveries += 1
        vaccinesDelivered += deliveryQty
        pod.vaccinesInStock += deliveryQty
        thisDelivery = vaccineDelivery(deliveryQty * 1, t + mdP)
        pod.vaccineDeliveries.append(thisDelivery)
        ratios[podIndex] = flightPreventedExposures(pod, droneVC, params, 1) / pod.flightTime
        displayTime = str(int((droneAvail[droneIndex])/60 + 7)) + ":" + "{:02d}".format((int(droneAvail[droneIndex]))%60)
        #print(displayTime,"- Drone", droneIndex+1, " delivered", deliveryQty,"vaccines to", pod.name)
        time += 5
    return deliveries, vaccinesDelivered
   
def findAvailDrone(droneAvail, time, pod, workMins):
    mI = np.argmin(droneAvail)
    if droneAvail[mI] <= time and time + pod.flightTime <= workMins:
        return mI
    else:
        return -1
   
def deliverByDrone(strategy, t, PODs, workingMinutes, droneVC, numDrones, params, mdP):
    if strategy == 'EPE':     #clever exposure-prevention delivery schedule
        delivDetails = deliverVaccinesByDroneEPE(t, PODs, workingMinutes, droneVC, numDrones, params, mdP)
    else:                   #simple delivery schedule
        delivDetails = deliverVaccinesByDroneSimple(strategy, t, PODs, workingMinutes, droneVC, numDrones, params, mdP)
    return delivDetails
      
def deliverByVehicle(openRoads, roadTimes, DCindex, strategy, PODs, workMinsPerDay, 
                     vehicleCapacities, numVehicles, params, mdP, t, maxTripLength, 
                     deliveryType):
    #Use Clarke-Wright method to find vehicle routes that could be used
    routes, routeTimes = findVehicleRoutes(openRoads, roadTimes, DCindex)
   
    displayTime = "7:00"
    totalVDelivered = 0
    times = np.zeros(numVehicles)
    
    vehiclesFinished = np.full(numVehicles, False)
    while all(vehiclesFinished) == False:
        for vehicle in range(0,numVehicles):
            if vehiclesFinished[vehicle] == True:
                continue
            #calculate the payoff for each route, acc. to the strategy
            routePayoffs = []
            for i in range(0,len(routes)):
                if times[vehicle] >= workMinsPerDay:
                    vehiclesFinished[vehicle] = True
                    routePayoffs.append(-1)
                    break
                #disqualify routes where less than 50 vaccines are needed per location, on average
                if calcVaccinesNeededOnRoute(routes[i], PODs, mdP) <= 50*(len(routes[i])-2):
                    routePayoffs.append(-1)
                    continue
                #disqualify routes (leaving after the first hour) that would exceed the day's working hours.
                if (times[vehicle] + routeTimes[i]) > workMinsPerDay and times[vehicle] > 60:
                    routePayoffs.append(-1)
                    continue
                #disqualify routes that take too long
                if routeTimes[i] > maxTripLength and deliveryType == "combined":
                    #the route is too long for the combined method
                    routePayoffs.append(-1)
                    continue
                
                if strategy == 'N':
                    totalN = 0
                    for node in routes[i]:
                        if node != DCindex:
                            totalN += PODs[node].N
                    routePayoffs.append(totalN / routeTimes[i])
                elif strategy == 'I':
                    totalI = 0
                    for node in routes[i]:
                        if node != DCindex:
                            totalI += PODs[node].I
                    routePayoffs.append(totalI / routeTimes[i])
                elif strategy == 'S':
                    totalS = 0
                    for node in routes[i]:
                        if node != DCindex:
                            totalS += PODs[node].S
                    routePayoffs.append(totalS / routeTimes[i])
                elif strategy == 'EPE':
                    #clever method, using expected exposures prevented
                    vehicleCapacity = vehicleCapacities[0]
                    if isRouteRural(routes[i], PODs):
                        vehicleCapacity = vehicleCapacities[1]
                    
                    totalEPE = 0
                    for node in routes[i]:
                        if node != DCindex:
                            vaccsNeeded = calcMaxVaccinesNeeded50(PODs[node], mdP)
                            vaccs = min(vaccsNeeded, vehicleCapacity)                                        
                            totalEPE += vehiclePreventedExposures(PODs[node], vaccs, params)
                    routePayoffs.append(totalEPE / routeTimes[i])
                else:
                    print("Invalid selection of strategy:", strategy)
                
            routeTimes = np.array(routeTimes)
            routePayoffs = np.array(routePayoffs)
            
            if max(routePayoffs) == -1:
                vehiclesFinished[vehicle] = True
                break
                        
            chosenRouteI = np.argmax(routePayoffs)
            
            #setting vehicle capacity for delivery
            vehicleCapacity = vehicleCapacities[0]  #deliver by default in large carrier
            if isRouteRural(routes[chosenRouteI], PODs):
                #deliver in a handheld carrier
                vehicleCapacity = vehicleCapacities[1]
                    
            vDelivered, strippedRoute = deliverRoute(routes[chosenRouteI], PODs, strategy, vehicleCapacity, mdP, t, params)
            
            totalVDelivered += sum(vDelivered)
            times[vehicle] += routeTimes[chosenRouteI]
            routeNames = []
            for location in strippedRoute:
                routeNames.append(PODs[location].name)
            
            displayTime = str(int(times[vehicle]/60 + 7)) + ":" + "{:02d}".format(int(times[vehicle]%60))
            #print(displayTime, "- Vehicle", vehicle+1, "delivered quantities", vDelivered, "to route", routeNames)
            
    #print(displayTime,"Day deliveries complete.")
    totMinsDriven = sum(times)
    return totalVDelivered, totMinsDriven

def isRouteRural(route, PODs):
    anyRurals = False
    for node in route:
        if PODs[node].isUrban == False:
            anyRurals = True
    return anyRurals
                
def deliverRoute(totalRoute, PODs, strategy, vehicleCapacity, mdP, t, params):
    '''Delivers vaccines to the PODs in the route, according to the amount still required at each POD, 
    with <strategy> determining which POD gets allocated the most vaccines.'''
    route = totalRoute[1:-1]        #removes DC from route string
    vaccinesToDeliver = np.zeros(np.shape(route))
    routeSINC = np.zeros(np.shape(route))   #this routeSINC is the S, I, N or C value for each 
    routeNeeds = np.zeros(np.shape(route))  #this routeNeeds is the number of vaccines each location needs
    for i in range(0,len(route)):
        routeNeeds[i] = calcMaxVaccinesNeeded50(PODs[route[i]], mdP)
        if strategy == 'S':
            routeSINC[i] = PODs[route[i]].S
        elif strategy == 'I':
            routeSINC[i] = PODs[route[i]].I
        elif strategy == 'N':
            routeSINC[i] = PODs[route[i]].N
        elif strategy == 'EPE':
            #for each pod in the route, calculate the prevented exposures if its given the max vaccines it needs
            routeSINC[i] = vehiclePreventedExposures(PODs[route[i]], routeNeeds[i], params)
    
    if sum(routeNeeds) < vehicleCapacity:
        vaccinesToDeliver = routeNeeds
    else:
        remainingCap = vehicleCapacity
        for i in range(0,len(route)):
            bestSINCindex = np.argmax(routeSINC)
            vaccinesToDeliver[bestSINCindex] = min(routeNeeds[bestSINCindex], remainingCap)
            routeSINC[bestSINCindex] = -1                       #do not add more vaccines to this spot.
            remainingCap -= vaccinesToDeliver[bestSINCindex]
            if remainingCap == 0:
                break
        
    for i in range(0,len(route)):
        vaccQty = vaccinesToDeliver[i]
        if vaccQty <= 0:
            continue
        podIndex = route[i]
        PODs[podIndex].vaccinesInStock += vaccQty
        thisDelivery = vaccineDelivery(vaccQty * 1, t + mdP)
        PODs[podIndex].vaccineDeliveries.append(thisDelivery)
    return vaccinesToDeliver, route
        
def calcVaccinesNeededOnRoute(totalRoute, PODs, mdP):
    ''' Given the route, PODs array and number of days the vaccine is potent for, 
    returns the expected number of vaccines needed to be delivered on the route.'''
    route = totalRoute[1:-1]
    totVaccsNeeded = 0
    for podIndex in route:
        totVaccsNeeded += calcMaxVaccinesNeeded50(PODs[podIndex], mdP)
    return totVaccsNeeded
       
def choosePODtoFlyTo(strategy, PODs, time, workingMinutesPerDay, mdP):
    '''Selects the POD with the highest number of susceptible people, among 
    PODs with less than 2000 vaccines. Returns -1 if all PODs have enough vaccines.'''
    maxV = -1
    maxpod = -1
    for pod in PODs:
        #Different strategies for selecting the order of deliveries
        if pod.flightTime == 0:
            continue
        elif strategy == 'S':
            podVal = pod.S / pod.flightTime
        elif strategy == 'I':
            podVal = pod.I / pod.flightTime
        elif strategy == 'N':
            podVal = pod.N / pod.flightTime
        else:
            print("Invalid delivery strategy selected. Choose S, I, N, or EPE.")
            quit()
        
        #the drone flight must be able to return in time
        if time + pod.flightTime <= workingMinutesPerDay:
            #there's no need for more vaccines than people to vaccinate
            if calcMaxVaccinesNeeded50(pod, mdP) > 0:
                #no drone flights to DC, and vaccine stock can't exceed max stock
                if podVal > maxV and pod.flightTime != 0:
                    maxV = podVal
                    maxpod = pod
    
    if maxV == 0:
        return -1
    return maxpod
                
def vehiclePreventedExposures(pod, vaccQty, params):
    #create two variations of this POD, with differing vaccinesInStock values
    beforePOD = copy.deepcopy(pod)
    afterPOD = copy.deepcopy(pod)
    afterPOD.vaccinesInStock += vaccQty
    thisDelivery = vaccineDelivery(vaccQty, np.inf)
    afterPOD.vaccineDeliveries.append(thisDelivery)
    
    #vaccinate both PODs with the vaccines available at each
    vaccinateOnePOD(beforePOD)
    vaccinateOnePOD(afterPOD)
    
    #progress both PODs another day
    progressSinglePOD(beforePOD, params)
    progressSinglePOD(afterPOD, params)
    
    #return the difference between exposure values - how many exposures did the flights prevent
    return max(beforePOD.E - afterPOD.E, 0)
                
def flightPreventedExposures(pod, droneVC, params, numFlights):
    #create two variations of this POD, with differing vaccinesInStock values
    beforeFlightsPOD = copy.deepcopy(pod)
    afterFlightsPOD = copy.deepcopy(pod)
    afterFlightsPOD.vaccinesInStock += droneVC * numFlights
    thisDelivery = vaccineDelivery(droneVC * numFlights, np.inf)
    afterFlightsPOD.vaccineDeliveries.append(thisDelivery)
    
    #vaccinate both PODs with the vaccines available at each
    vaccinateOnePOD(beforeFlightsPOD)
    vaccinateOnePOD(afterFlightsPOD)
    
    #progress both PODs another day
    progressSinglePOD(beforeFlightsPOD, params)
    progressSinglePOD(afterFlightsPOD, params)
    
    #return the difference between exposure values - how many exposures did the flights prevent
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
    
def findVehicleRoutes(oR, rT, dcI):
    '''Finds routes of minimum time taken, using Clarke-Wright algorithm.
    Also includes routes to-and-from each POD using Dijkstra's shortest path algorithm.'''
    openRoadTimes = np.full(np.shape(oR), 100000)   #BigM if no direct road from i to j. Otherwise, travel time in mins
    for i in range(0,len(openRoadTimes)):
        for j in range(0,len(openRoadTimes[i])):
            if oR[i][j] == 1 and rT[i][j] != np.inf:
                openRoadTimes[i][j] = rT[i][j]
                if openRoadTimes[i][j] == np.inf:
                    openRoadTimes[i][j] = 100000    #BigM easier to work with than infinity
    #Calculate the savings travelling directly from each pair to each other produces
    Savings = calcSavings(openRoadTimes, dcI)
    #build up the vehicle routes according to the heuristic
    routes = buildRoutes(Savings, dcI, openRoadTimes)
    routeTimes = []
    for route in routes:
        routeTimes.append(calcRouteTime(route, openRoadTimes))
        
    #Add to the routes found, there-and-back routes for each possible POD (back to DC).
    for i in range(0, len(openRoadTimes)):
        if i == dcI:
            continue
        route = [dcI,i,dcI]
        rTime = 2 * findShortestPath(dcI, i, openRoadTimes)[0]
        if rTime < 100000:          #if the route is feasible, add it
            routes.append(route)
            routeTimes.append(rTime)   
    return routes, routeTimes
    
def calcSavings(oRT, dcI):
    S = np.zeros(np.shape(oRT))
    for i in range(0,len(oRT)):
        for j in range(0,len(oRT[i])):
            if i == dcI or j == dcI:
                S[i][j] = -100000
                continue
            if i != j:
                S[i][j] = oRT[dcI][i] + oRT[dcI][j] - oRT[i][j]        
    return S            
    
def buildRoutes(Savings, dcI, openRoadTimes):
    routes = []    
    inRoutes = np.full(len(Savings), False)
    while np.amax(Savings) > 0:
        #Select savings pairs in descending order to add to routes
        rowMaxes = np.amax(Savings, axis=1)
        #rowI, colI is the location pair with the highest savings
        rowI = np.argmax(rowMaxes)
        colI = np.where(rowMaxes[rowI] == Savings[rowI])[0][0]
        Savings[rowI][colI] = 0
        Savings[colI][rowI] = 0
        if openRoadTimes[rowI][colI] > 90000:
            #only allow this link if feasible (capacity), etc
            continue
        
        if inRoutes[rowI] == inRoutes[colI] == False:
            routes.append([rowI, colI])
            inRoutes[rowI] = inRoutes[colI] = True
        elif inRoutes[rowI] == inRoutes[colI] == True:
            RrouteIndex = findPoint(routes, rowI, dcI)
            CrouteIndex = findPoint(routes, colI, dcI)
            #Neither POD is interior, and the PODs are not in the same route
            if RrouteIndex != -1 and CrouteIndex != -1 and CrouteIndex != RrouteIndex:
                if routes[RrouteIndex][0] == rowI:      #rowI must be last in its route
                    routes[RrouteIndex].reverse()
                if routes[CrouteIndex][-1] == colI:     #colI must be first in its route
                    routes[CrouteIndex].reverse()
            
                routes[RrouteIndex].extend(routes[CrouteIndex])
                if RrouteIndex != CrouteIndex:  
                    del routes[CrouteIndex]
        else:
            if inRoutes[rowI] == True:      #the rowI is in a route already, colI not
                routeIndex = findPoint(routes, rowI, dcI)
                if routeIndex == -1:
                    continue
                
                if routes[routeIndex][0] == rowI:
                    routes[routeIndex].reverse()
                    routes[routeIndex].append(colI)
                elif routes[routeIndex][-1] == rowI:
                    routes[routeIndex].append(colI)
                else:
                    print("Logical error.")
                    quit()
                inRoutes[colI] = True
            else:                           #the colI is in a route already, rowI not
                routeIndex = findPoint(routes, colI, dcI)
                if routeIndex == -1:
                    continue
                
                if routes[routeIndex][0] == colI:
                    routes[routeIndex].reverse()
                    routes[routeIndex].append(rowI)
                elif routes[routeIndex][-1] == colI:
                    routes[routeIndex].append(rowI)
                else:
                    quit()
                inRoutes[rowI] = True
                
    #Need to add the visit to the DC at the start and end of each
    for i in range(0,len(routes)):
        routes[i].append(dcI)
        routes[i].reverse()
        routes[i].append(dcI)
        routes[i].reverse()
        
    return routes
    
def findPoint(routes, node, dcI):
    '''If the node is interior or the DC, returns -1. If the node is on an edge,
    returns the index of the route it is in.'''
    if node == dcI:
        return -1
    
    routeIndex = -1
    for route in routes:
        routeIndex += 1
        if route[0] == node or route[-1] == node:
            return routeIndex
    return -1
  
def calcRouteTime(route, oRT):
    '''Calculates the total time in minutes to travel the given route.
    If the route contains any indirect links between nodes (where there is no passable road), 
    use Dijkstra's shortest path algorithm to determine to time to traverse that link if possible.'''
    routeTime = 0
    for k in range(1, len(route)):
        i = route[k-1]
        j = route[k]
        if oRT[i][j] == 100000:       #indirect link between nodes i and j
            ijTime, ijPath = findShortestPath(i,j,oRT)
            #ijTime is the time of the shortest path from i to j. ijPath, the actual shortest path.
            routeTime += ijTime
        else:                           #add time to travel direct link
            routeTime += oRT[i][j]
    return routeTime
  
def findShortestPath(i,j,oRT):
    '''Finds the shortest path between i and j using Dijkstra's algorithm.'''  
    nodeVisited = np.full(len(oRT), False)  #True if node has been visited in Dijkstra's alg
    times = np.full(len(oRT), 100000)       #holds the min time to reach each node from initial
    prevNode = np.full(len(oRT), -1)        #holds the previous node in the shortest path
    times[i] = 0
    while all(nodeVisited) == False:
        minUnvisitedTime = 100000
        minUnvisitedNode = -1
        for i in range(0,len(nodeVisited)):
            if nodeVisited[i] == False:
                if times[i] < minUnvisitedTime:
                    minUnvisitedTime = times[i]
                    minUnvisitedNode = i
        if minUnvisitedNode == -1:
            break
        currNode = minUnvisitedNode
        nodeVisited[currNode] = True
        #update distance values
        for node in range(0,len(oRT[currNode])):
            if oRT[currNode][node] != 100000 and oRT[currNode][node] != 0:
                if times[currNode] + oRT[currNode][node] < times[node]:
                    times[node] = times[currNode] + oRT[currNode][node]
                    prevNode[node] = currNode
    pN = -2
    ijPath = [j]
    while pN != -1:
        pN = prevNode[ijPath[-1]]
        if pN != -1:
            ijPath.append(pN)
    ijPath.reverse()
    return times[j], ijPath
  
def assignTeams(PODs, numTeams, tStrategy):
    if tStrategy == 'spread':
        remainingTeams = numTeams
        for pod in PODs:
            pod.teamsAtPOD = min(roundUsingProb(numTeams/len(PODs)), remainingTeams)
            remainingTeams -= pod.teamsAtPOD
            pod.maxVaccinationsPerDay = pod.vaccsPerTeamDay * pod.teamsAtPOD
            
        if remainingTeams > 0:
            randIndex = int(len(PODs) * random.random())
            PODs[randIndex].teamsAtPOD += remainingTeams
            remainingTeams = 0
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
        
    plt.plot(np.arange(daysOfIntervention), Ss, label='Susceptible')
    plt.plot(np.arange(daysOfIntervention), Es, label='Exposed')
    plt.plot(np.arange(daysOfIntervention), Is, label='Infectious')
    plt.plot(np.arange(daysOfIntervention), deads, label='Deaths')
    plt.plot(np.arange(daysOfIntervention), Rs, label='Recovered')
    plt.plot(np.arange(daysOfIntervention), (vacs-min(vacs)), label='Vaccinations')

    label = ""
    for i in podIndexes:
        label = label + PODs[i].name + ","
    label = label[:len(label)-1]
    if len(podIndexes) == len(PODs):
        label = "Measles epidemic progression over the whole network,\nwith intervention"
    
    plt.title(label)
    #plt.title("Progression of measles epidemic SEIR model")
    plt.ylabel("Number of People")
    plt.xlabel("Time elapsed in days")
    plt.legend()
    
#     import matplotlib as mpl
#     mpl.rcParams['figure.dpi'] = 150
#     plt.savefig("SEIRvehicleD.pdf",bbox_inches='tight')
    
    plt.show()
    
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
    PODs, roadDistances = readPODsFromFile(filename, maxVaccsFreePODday, maxVaccsFixedPODday, ruralTnt, urbanTnt)
    podDistances = calcPODdistanceMatrix(PODs)
    numRoadsClosed, openRoads, roadTimes = calcRoadTimes(PODs, roadDistances, podDistances, roadCloseFactor)
    #roadTimes = roadDistances   #for matadi only!
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
    totMinsDriven = 0
    totVaccsGiven = 0
    vaccsPerDay = np.zeros(simulationRuntime)

    #===============================================================================
    # Simulation - Main Loop
    #===============================================================================    
    for t in range(0, simulationRuntime):
        #print("\n",t,":")
            
        #Progress the SEIR models and migration by one day
        PODs = progressEpidemicByOneDay(PODs, params, MigrationProportions)
        
        #if the intervention is running today
        if t >= interventionStartTime and t < interventionStartTime + interventionLength and random.random() < workDaysPerWeek/7:
            waitingForIntervention = False
            
            #expire old vaccines
            totExpired += expireVaccines(PODs, t)
            
            #assign team locations
            assignTeams(PODs, numTeams, teamStrategy)
            
            #deliver vaccines
            if deliveryType == "drone":
                delivDetails = deliverByDrone(strategy, t, PODs, workingMinutesPerDay, 
                                            droneVaccineCapacity, numberOfDrones, params, 
                                            monoDaysPotency)
                totDroneDelvs += delivDetails[0]
                totVaccs += delivDetails[1]
            elif deliveryType == "vehicle":
                vDeliveredToday, minsDr = deliverByVehicle(openRoads, roadTimes, DCindex, 
                                            strategy, PODs, workingMinutesPerDay, 
                                            vehicleCapacities, numVehicles, params, 
                                            monoDaysPotency, t, maxTripLength, deliveryType)
                totVaccs += vDeliveredToday
                totMinsDriven += minsDr
            elif deliveryType == "combined":
                vaccsByV, minsDriven = deliverByVehicle(openRoads, roadTimes, DCindex, 
                                            strategy, PODs, workingMinutesPerDay, 
                                            vehicleCapacities, numVehicles, params, 
                                            monoDaysPotency, t, maxTripLength, deliveryType)
                dDelivs, vaccsByD = deliverByDrone(strategy, t, PODs, workingMinutesPerDay, 
                                            droneVaccineCapacity, numberOfDrones, params, 
                                            monoDaysPotency)
                totVaccs += vaccsByD + vaccsByV
                totDroneDelvs += dDelivs
                totMinsDriven += minsDriven
            elif deliveryType != "none":            #no vaccinations happen if deliveryType is none.
                print("Invalid delivery type selected.")
                quit()        
            
            #process vaccinations for the day
            if deliveryType != "none":
                vGiven, vGivenDC = vaccinate(PODs)
                vaccsPerDay[t] = vGiven
                totVaccsGiven += vGiven
                totVaccs += vGivenDC
        elif t >= interventionStartTime + interventionLength:
            #intervention is over
            interventionOver = True
            #print("Intervention is over")
        elif not waitingForIntervention and t <= interventionStartTime:
            for pod in PODs:
                if pod.isUrban and pod.I/pod.N > interventionCaseRatio:
                    #if a confirmed measles case found in a major town, then decide to intervene.
                    interventionStartTime = t + interventionLeadTime
                    waitingForIntervention = True
                    #print("The epidemic has been detected in", pod.name, "! Intervention will begin on day", interventionStartTime)            
                    break
        
        drivingCost = totMinsDriven/60 * 60 * 7/100 * 1.36  #60kph, 7l/100km, $1.36 per litre
        deliveryCost = totDroneDelvs * costPerFlight + drivingCost
        vaccineCost = totVaccs * costPerDoseMono
            
        updatePlotHistory(PODs, plots)
        
        #plot this iteration's png
        #plotMap(PODs, t, waitingForIntervention, interventionStartTime, interventionOver, 
        #       totExpired, totVaccs, totVaccsGiven, deliveryCost + vaccineCost)
    
    #===============================================================================
    # Result Reporting
    #===============================================================================
    deaths = 0
    for pod in PODs:
        deaths += pod.deaths
    #print("\nThe total number of deaths is:", deaths)
    
    if deliveryType == "drone":
        droneCost = totDroneDelvs * costPerFlight
        #print("Total drone flight cost: $", droneCost,",for", totDroneDelvs, "flights.")
    elif deliveryType == "vehicle":
        drivingCost = totMinsDriven/60 * 60 * 7/100 * 1.36     #60kph, 7l/100km, $1.36 per litre
        #print("Total estimated driving cost: $", int(drivingCost), "for a total of", int(totMinsDriven/60), "hours of driving.")    
    elif deliveryType == "combined":
        drivingCost = totMinsDriven/60 * 60 * 7/100 * 1.36     #60kph, 7l/100km, $1.36 per litre
        #print("Total estimated driving cost: $", "for a total of", int(totMinsDriven/60), "hours of driving.")    
        droneCost = totDroneDelvs * costPerFlight
        #print("Total drone flight cost: $", droneCost,",for", totDroneDelvs, "flights.")
    
    
    #print("Total cost of monodose vaccines delivered: $", vaccineCost, ",for", totVaccs, "doses.")
    if totVaccs > 0:
        expiryRatio = totExpired / totVaccs * 100
        #print("Percentage of vaccines expired without use:", round(expiryRatio,2), "%")
    
    #print("Total number of vaccines administered to patients:", totVaccsGiven)
    
    #plotPODSum(simulationRuntime, plots, (0,1,2,3,4,5,6,7,8,9,10), PODs)
    #plotPODSum(simulationRuntime, plots, (0,1,2,3,4,5,6,7,8,9,10,11), PODs)
    #plotPOD(simulationRuntime, plots, 7, "Likasi")
    #plotMap(PODs,  t, waitingForIntervention, interventionStartTime, interventionOver,
     #      totExpired, totVaccs, totVaccsGiven, vaccineCost + deliveryCost)
    
    Vactots = np.zeros(simulationRuntime)
    Stots = np.zeros(simulationRuntime)
    Itots = np.zeros(simulationRuntime)
    for t in range(0,simulationRuntime):    
        Vactots[t] += plots[6][5][t]  #Just for Likasi
        Stots[t] += plots[6][0][t]  
        Itots[t] += plots[6][2][t]
#         for podI in range(0,len(PODs)):        #for all locations
#             Vactots[t] += plots[podI][5][t]  
#             Stots[t] += plots[podI][0][t]  
#             Itots[t] += plots[podI][2][t]
    return deaths, vaccineCost + deliveryCost, totVaccsGiven, totExpired#, vaccsPerDay #,Vactots,Stots


simulationRuntime = 150             #days to run the simulation for
#Parameters    ========================================================================
#Measles SEIR parameters
exposedDays = 10                    #number of days a patient is exposed for without symptoms
infectiousDays = 8                  #number of days a patient is infectious for
deathRate = 0.10 * 1/infectiousDays #the proportion of infected patients that die per day
R0 = 15                             #basic reproductive number of the epidemic
params = [R0 / infectiousDays, 1/exposedDays, 1/infectiousDays, deathRate]
migrationIntensity = 1              #factor by which migration is multiplied. 2 means more migration.
#vaccine parameters
vaccineEffectiveness = 0.95         #probability the vaccine works (for non-exposed)
prophylaxis72hrSuccessRate = 0.83   #probability the vaccine works (for exposed, within 72hrs)
monoDaysPotency = 3                 #number of days for which the vaccine lasts outside of cold-chain
TenDaysPotency = 5                  #number of days the 10-dose vaccine lasts outside of cold-chain
#intervention parameters
interventionLeadTime = 15           #number of days before vaccination starts
interventionCaseRatio = 0.005       #ratio of I/S in a town before detection
interventionLength = 28             #number of days the intervention lasts for
workingMinutesPerDay = 660          #11 working hours per day: 7am to 6pm
workDaysPerWeek = 7                 #number of working days per week for MSF teams
numTeams = 15                       #number of vaccination teams in the field
maxVaccsFixedPODday = 450           #fixed teams can vaccinate up to 450 per day
maxVaccsFreePODday = 250            #free teams can vaccinate up to 250 per day
ruralTnt = 450                      #rural areas, turnout is 300-600
urbanTnt = 900                      #urban areas, turnout is 800-1000
#delivery details
flightLaunchTime = 10               #minutes per flight, to set up takeoff
droneSpeed = 100                    #100 kilometres per hour     
numberOfDrones = 2                  #number of drones
numVehicles = 2                     #number of land-based delivery vehicles
droneVaccineCapacity = 60           #number of vaccine doses per drone
vehicleCapacities = [1050,200]      #number of vaccine doses per vehicle
roadCloseFactor = 1                 #prob that road is open is divided by this factor. 1<rCF
#costs
costPerDoseMono = 2.85              #the cost per dose of monodose measles vaccine
costPerDose10 = 1.284               #the cost per dose of 10-dose measles vaccine
costPerFlight = 17                  #$17 per drone flight
#strategies
strategy = 'I'                      #I = infections, S = Susceptible, N = Total Pop., EPE = Expected Prevented Exposures
teamStrategy = 'N'                  #I, S, N, I/N, spread
deliveryType = 'vehicle'            #"none", "drone" or "vehicle", or "combined"
maxTripLength = 180                 #if deliveryType = combined, this is the cutoff in mins for vehicle trip length

#changes made in calcMaxVaccinesNeeded50, and both deliverByDrones methods, using np.ceil to round up to multips of 60

print(simulate())