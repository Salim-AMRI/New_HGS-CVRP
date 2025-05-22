void Genetic::run()
{
        /* INITIAL POPULATION */
        population.generatePopulation();

        int cpt = 1;

        std::ofstream Evolutionfile;
        if (params.ap.crossoverType == 10) {
            Evolutionfile.open("evolution/Instance_" + std::to_string(params.nbClients + 1) + "_" + "crossover_type_" + std::to_string(params.ap.crossoverType) + "_nbCut_" + std::to_string(params.ap.nbCut) + "_" + "eqSeg_" + std::to_string(params.ap.eqSeg) + "_" + "useCostBenefit_" + std::to_string(params.ap.useCostBenefit) + "_" + "randSelect_" + std::to_string(params.ap.randSelect) + "_" + "insertSeg_" + std::to_string(params.ap.insertSeg) + "_seed_" +  std::to_string(params.ap.seed) + ".txt");
        } else {
            Evolutionfile.open("evolution/Instance_" + std::to_string(params.nbClients + 1) + "_" + "crossover_type_" + std::to_string(params.ap.crossoverType) + "_seed_" +  std::to_string(params.ap.seed) + ".txt");
        }
        
        const double interval = 0.2; // Interval de temps en secondes
        double lastOutputTime = (double)clock() / (double)CLOCKS_PER_SEC; // Initialisation du temps de la dernière sortie

        int nbIter;
        int nbIterNonProd = 1;
        if (params.verbose) std::cout << "----- STARTING GENETIC ALGORITHM" << std::endl;

        std::cout << " params.ap.nbIter " <<   params.ap.nbIter << std::endl;
        std::cout << " params.ap.crossoverType " <<   params.ap.crossoverType << std::endl;

        for (nbIter = 0 ;  (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC < params.ap.timeLimit ; nbIter++)
        {
                /* SELECTION AND CROSSOVER */
                if (params.ap.crossoverType == 0) {
                    crossoverOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 1) {
                    crossoverLOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 2) {
                    crossoverAOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 3) {
                    crossoverUOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 4) {
                    crossoverPMX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 5) {
                    crossoverEAX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 6) {
                    crossoverGPX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 7) {
                    crossoverGrPX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 8) {
                    crossoverPathRelinking(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 9) {
                    crossoverGrPX2(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 10) {
                    crossoverGOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                }

        
        
        double currentTime = (double)clock() / (double)CLOCKS_PER_SEC;
        if (currentTime - lastOutputTime >= interval)
        {
            Evolutionfile << population.getBestFeasible()->eval.penalizedCost << " - " << currentTime << std::endl;
            lastOutputTime = currentTime;
        }


        /* LOCAL SEARCH */
        localSearch.run(offspring, params.penaltyCapacity, params.penaltyDuration);
                bool isNewBest = population.addIndividual(offspring,true);
                if (!offspring.eval.isFeasible && params.ran()%2 == 0) // Repair half of the solutions in case of infeasibility
                {
                        localSearch.run(offspring, params.penaltyCapacity*10., params.penaltyDuration*10.);
                        if (offspring.eval.isFeasible) isNewBest = (population.addIndividual(offspring,false) || isNewBest);
                }

                /* TRACKING THE NUMBER OF ITERATIONS SINCE LAST SOLUTION IMPROVEMENT */
                if (isNewBest) nbIterNonProd = 1;
                else nbIterNonProd ++ ;

                /* DIVERSIFICATION, PENALTY MANAGEMENT AND TRACES */
                if (nbIter % params.ap.nbIterPenaltyManagement == 0) population.managePenalties();
                if (nbIter % params.ap.nbIterTraces == 0) population.printState(nbIter, nbIterNonProd);

                /* FOR TESTS INVOLVING SUCCESSIVE RUNS UNTIL A TIME LIMIT: WE RESET THE ALGORITHM/POPULATION EACH TIME maxIterNonProd IS ATTAINED*/

// 		if (params.ap.timeLimit != 0 && nbIterNonProd == params.ap.nbIter)
// 		{
// 			population.restart();
// 			nbIterNonProd = 1;
//             MyFile << "Restart" << std::endl;
//
// 		}

        }


        Evolutionfile << population.getBestFeasible()->eval.penalizedCost << " - " << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << std::endl;

        if (params.verbose) std::cout << "----- GENETIC ALGORITHM FINISHED AFTER " << nbIter << " ITERATIONS. TIME SPENT: " << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << std::endl;

        Evolutionfile.close();
}

void Genetic::crossoverOX(Individual & result, const Individual & parent1, const Individual & parent2)
{
        // Frequency table to track the customers which have been already inserted
        std::vector <bool> freqClient = std::vector <bool> (params.nbClients + 1, false);

        // Picking the beginning and end of the crossover zone
        std::uniform_int_distribution<> distr(0, params.nbClients-1);
        int start = distr(params.ran);
        int end = distr(params.ran);

        // Avoid that start and end coincide by accident
        while (end == start) end = distr(params.ran);

        // Copy from start to end
        int j = start;
        while (j % params.nbClients != (end + 1) % params.nbClients)
        {
                result.chromT[j % params.nbClients] = parent1.chromT[j % params.nbClients];
                freqClient[result.chromT[j % params.nbClients]] = true;
                j++;
        }

        // Fill the remaining elements in the order given by the second parent
        for (int i = 1; i <= params.nbClients; i++)
        {
                int temp = parent2.chromT[(end + i) % params.nbClients];
                if (freqClient[temp] == false)
                {
                        result.chromT[j % params.nbClients] = temp;
                        j++;
                }
        }

        // Complete the individual with the Split algorithm
        split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverLOX(Individual & result, const Individual & parent1, const Individual & parent2)
{
    // Frequency table to track the customers which have been already inserted
    std::vector<bool> freqClient = std::vector<bool>(params.nbClients + 1, false);

    // Picking the beginning and end of the crossover zone
    std::uniform_int_distribution<> distr(0, params.nbClients - 1);
    int start = distr(params.ran);
    int end = distr(params.ran);

    // Ensure start is smaller than end
    if (start > end)
        std::swap(start, end);

    // Copy the crossover segment from parent1 to the offspring
    for (int i = start; i <= end; i++) {
        result.chromT[i % params.nbClients] = parent1.chromT[i % params.nbClients];
        freqClient[result.chromT[i % params.nbClients]] = true;
    }

    // Fill the remaining elements using parent2
    int j = (end + 1) % params.nbClients; // Start from the next index after end
    for (int i = 0; i < params.nbClients; i++) {
        if (!freqClient[parent2.chromT[i]]) {
            result.chromT[j] = parent2.chromT[i];
            freqClient[parent2.chromT[i]] = true;
            j = (j + 1) % params.nbClients;
        }
    }

    // Complete the individual with the Split algorithm
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverAOX(Individual & result, const Individual & p1, const Individual & p2)
{
    // Frequency table to track customers already inserted
    std::vector<bool> freqClient = std::vector<bool>(params.nbClients + 1, false);

    // Picking the beginning and end of the crossover zone
    int start1 = std::rand() % params.nbClients;
    int end1 = std::rand() % params.nbClients;
    while (end1 == start1) end1 = std::rand() % params.nbClients;

    // Shift zone in p2 to match final customer of zone in p1
    int start2 = start1, end2 = end1;
    while (p2.chromT[end2 % params.nbClients] != p1.chromT[end1 % params.nbClients])
        start2++, end2++;

    // Test if zone in p1 is different to zone in p2
    bool same = true;
    int size = (start1 < end1 ? end1 - start1 : params.nbClients - start1 + end1);
    for (int j = 0; j < size && same; j++) {
        if (p1.chromT[(start1 + j) % params.nbClients] != p2.chromT[(start2 + j) % params.nbClients])
            same = false;
    }

    // If same, randomize point in p2
    if (same) end2 = end2 + rand() % (params.nbClients - size);

    // Copy in place the elements from start to end
    int j = start1;
    while (j % params.nbClients != (end1 + 1) % params.nbClients) {
        result.chromT[j % params.nbClients] = p1.chromT[j % params.nbClients];
        freqClient[result.chromT[j % params.nbClients]] = true;
        j++;
    }

    // Fill the remaining elements in the order given by p2
    for (int i = 1; i <= params.nbClients; i++) {
        int temp = p2.chromT[(end2 + i) % params.nbClients];
        if (freqClient[temp] == false) {
            result.chromT[j % params.nbClients] = temp;
            j++;
        }
    }

    // Completing the individual with the Split algorithm
    split.generalSplit(result, p1.eval.nbRoutes);
}

void Genetic::crossoverUOX(Individual & result, const Individual & parent1, const Individual & parent2)
{
    // Frequency table to track the customers which have been already inserted
    std::vector<bool> freqClient = std::vector<bool>(params.nbClients + 1, false);

    // Picking crossover points
    std::uniform_int_distribution<> distr(0, params.nbClients - 1);
    int point1 = distr(params.ran);
    int point2 = distr(params.ran);

    // Ensure point1 and point2 are different
    while (point2 == point1)
        point2 = distr(params.ran);

    // Ensure point2 is greater than point1
    if (point1 > point2)
        std::swap(point1, point2);

    // Copy elements between points from parent1 to result
    for (int i = point1; i <= point2; ++i) {
        result.chromT[i] = parent1.chromT[i];
        freqClient[result.chromT[i]] = true;
    }

    // Fill remaining elements from parent2
    int j = (point2 + 1) % params.nbClients;
    for (int i = 0; i < params.nbClients; ++i) {
        int index = (point2 + 1 + i) % params.nbClients;
        int temp = parent2.chromT[index];
        if (!freqClient[temp]) {
            result.chromT[j] = temp;
            freqClient[temp] = true;
            j = (j + 1) % params.nbClients;
        }
    }

    // Complete the individual with the Split algorithm
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverPMX(Individual & result, const Individual & parent1, const Individual & parent2)
{
    // Picking the beginning and end of the crossover zone
    std::uniform_int_distribution<> distr(0, params.nbClients - 1);
    int start = distr(params.ran);
    int end = distr(params.ran);

    // Avoid that start and end coincide by accident
    while (end == start)
        end = distr(params.ran);

    // Copy the segment from parent1 to the result and create mapping table
    for (int i = start; i <= end; ++i) {
        result.chromT[i] = parent1.chromT[i];
    }

    // Create mapping table for elements in the segment
    std::unordered_map<int, int> mapping;
    for (int i = start; i <= end; ++i) {
        mapping[parent1.chromT[i]] = parent2.chromT[i];
    }

    // Map elements from parent2 to result, avoiding elements already in the segment
    for (int i = 0; i < params.nbClients; ++i) {
        if (i < start || i > end) {
            int value = parent2.chromT[i];
            while (mapping.find(value) != mapping.end()) {
                value = mapping[value];
            }
            result.chromT[i] = value;
        }
    }

    // Complete the individual with the Split algorithm
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverEAX(Individual & result, const Individual & parent1, const Individual & parent2)
{
    // Step 1: Create a AB consisting of all edges of parent1 and parent2
    std::vector<std::vector<int>> AB(params.nbClients, std::vector<int>(params.nbClients, 0));
    for (int i = 0; i < params.nbClients; ++i) {
        for (int j = 0; j < params.nbClients; ++j) {
            if (parent1.chromT[i] == j || parent2.chromT[i] == j) {
                AB[i][j] = 1;
            }
        }
    }

    // Step 2: Find closed cycles by tracing the edges of parent1 and parent2 alternately from GAB
    std::vector<bool> visited(params.nbClients, false);
    std::vector<std::vector<int>> cycles;
    for (int i = 0; i < params.nbClients; ++i) {
        if (!visited[i]) {
            std::vector<int> cycle;
            int currentNode = i;
            while (!visited[currentNode]) {
                cycle.push_back(currentNode);
                visited[currentNode] = true;
                int nextNode = -1;
                for (int j = 0; j < params.nbClients; ++j) {
                    if (!visited[j] && AB[currentNode][j] == 1) {
                        nextNode = j;
                        break;
                    }
                }
                if (nextNode == -1 || nextNode == i) {
                    break; // No unvisited neighbor found or cycle is complete
                }
                currentNode = nextNode;
            }
            if (!cycle.empty()) {
                cycles.push_back(cycle);
            }
        }
    }

    // Randomly shuffle the cycles
    std::random_shuffle(cycles.begin(), cycles.end());

    // Step 3: Construct Eset consisting of several AB-cycles selected randomly from all cycles
    std::vector<std::vector<int>> Eset;
    std::random_shuffle(cycles.begin(), cycles.end());
    for (const auto& cycle : cycles) {
        if (cycle.size() > 2) { // Exclure les cycles avec seulement deux arêtes
            // Vérifier et supprimer les clients ayant une valeur égale à 0
            std::vector<int> filteredCycle;
            for (int node : cycle) {
                if (node != 0) {
                    filteredCycle.push_back(node);
                }
            }
            if (!filteredCycle.empty()) {
                Eset.push_back(filteredCycle);
            }
        }
    }

    // Step 4: Manipulate edges of parent1 based on Eset
    for (const auto& cycle : Eset) {
        for (int i = 0; i < cycle.size(); ++i) {
            int currentNode = cycle[i];
            int nextNode = cycle[(i + 1) % cycle.size()];
            result.chromT[currentNode] = nextNode;
        }
    }

    // Step 5: Assemble sub-tours into a complete cycle
    bool completeCycle = false;
    while (!completeCycle) {
        int startNode = -1;
        for (int i = 0; i < params.nbClients; ++i) {
            if (!visited[i]) {
                startNode = i;
                break;
            }
        }
        if (startNode == -1) {
            completeCycle = true;
            break;
        }
        int currentNode = startNode;
        while (true) {
            visited[currentNode] = true;
            int nextNode = -1;
            for (int i = 0; i < params.nbClients; ++i) {
                if (!visited[i] && AB[currentNode][i] == 1) {
                    nextNode = i;
                    break;
                }
            }
            if (nextNode == -1) {
                result.chromT[currentNode] = startNode;
                break;
            }
            result.chromT[currentNode] = nextNode;
            currentNode = nextNode;
        }
    }

    // Déclaration d'un vecteur pour suivre les clients déjà vus dans le cycle final
    std::vector<bool> clientVus(params.nbClients, false);

    // Recherche des clients en double et leur remplacement
    for (int i = 0; i < params.nbClients; ++i) {
        int clientActuel = result.chromT[i];
        if (!clientVus[clientActuel]) {
            // Marquer le client comme vu
            clientVus[clientActuel] = true;
        } else {
            // Remplacer le client répété par le client manquant du parent1.chromT
            for (int j = 0; j < params.nbClients; ++j) {
                if (parent1.chromT[j] != 0 && !clientVus[parent1.chromT[j]]) {
                    result.chromT[i] = parent1.chromT[j];
                    clientVus[parent1.chromT[j]] = true;
                    break;
                }
            }
        }
    }

    // Complete the individual with the Split algorithm
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverGPX(Individual & result, const Individual & parent1, const Individual & parent2) {
    // Taille des chromosomes des parents
    int length = parent1.chromT.size();

    // Initialisation du chromosome de l'enfant résultant
    result.chromT.clear();

    // Définition de la matrice AB
    std::unordered_map<int, std::unordered_map<int, int>> AB;

    // Initialisation de la matrice AB
    for (int i = 0; i < length; ++i) {
        AB[i] = std::unordered_map<int, int>();
    }

    // Construction de la matrice AB qui combine les deux parents
    for (int i = 0; i < length - 1; ++i) {
        AB[parent1.chromT[i]][parent1.chromT[i + 1]] = 1; // Les arcs du parent 1
        AB[parent1.chromT[i + 1]][parent1.chromT[i]] = 1;

        // Vérifier si l'arête est partagée entre les deux parents
        int vertex1 = parent1.chromT[i];
        int vertex2 = parent1.chromT[i + 1];
        if (AB[parent2.chromT[i]][parent2.chromT[i + 1]] == 1 &&
            (parent2.chromT[i] == vertex1 || parent2.chromT[i] == vertex2) &&
            (parent2.chromT[i + 1] == vertex1 || parent2.chromT[i + 1] == vertex2)) {

            AB[vertex1][vertex2] = 2; // Les arcs partagés entre les deux parents
            AB[vertex2][vertex1] = 2;
        }
    }

    // Traitement pour le dernier arc du chromosome
    AB[parent1.chromT[length - 1]][parent1.chromT[0]] = 1;
    AB[parent1.chromT[0]][parent1.chromT[length - 1]] = 1;
    if (std::find(parent2.chromT.begin(), parent2.chromT.end(), parent1.chromT[length - 1]) != parent2.chromT.end() &&
        std::find(parent2.chromT.begin(), parent2.chromT.end(), parent1.chromT[0]) != parent2.chromT.end()) {
        AB[parent1.chromT[length - 1]][parent1.chromT[0]] = 2;
        AB[parent1.chromT[0]][parent1.chromT[length - 1]] = 2;
    }

    // Suppression des arêtes partagées (marquées avec la valeur 2) dans la matrice AB
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j) {
            if (AB[i][j] == 2) {
                AB[i][j] = -1; // Marquer temporairement pour la suppression
                AB[j][i] = -1; // Marquer temporairement pour la suppression
            }
        }
    }

    // Suppression des arêtes temporairement marquées
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j) {
            if (AB[i][j] == -1) {
                AB[i][j] = 0; // Suppression de l'arête
            }
        }
    }

    // Transformation de la matrice AB pour que tous les sommets aient un degré de 3
    for (int i = 0; i < length; ++i) {
        int degree = 0;
        for (int j = 0; j < length; ++j) {
            degree += AB[i][j]; // Calcul du degré du sommet i
        }
        for (int j = 0; j < length && degree < 3; ++j) {
            if (i != j && AB[i][j] == 0) {
                AB[i][j] = 1; // Ajouter une arête supplémentaire
                AB[j][i] = 1; // Arête bidirectionnelle
                degree++;
            }
        }
    }

    // Transformation de la matrice AB pour que tous les sommets aient un degré de 3
    int countDegreeThree = 0; // Initialisation du compteur
    for (int i = 0; i < length; ++i) {
        int degree = 0;
        for (int j = 0; j < length; ++j) {
            degree += AB[i][j]; // Calcul du degré du sommet i
        }
        if (degree == 3) {
            countDegreeThree++; // Incrémentation du compteur si le degré est égal à 3
        }
    }

    // Construction du chromosome de l'enfant résultant du crossover GPX
    std::vector<bool> visited(length, false); // Marque les sommets visités
    int current = 0; // Sommet de départ (peut être choisi arbitrairement)

    // Boucle de construction du chemin
    do {
        result.chromT.push_back(current); // Ajouter le sommet actuel à la séquence du chromosome de l'enfant
        visited[current] = true; // Marquer le sommet comme visité

        // Trouver le prochain sommet connecté non visité
        int next = -1;
        for (int i = 0; i < length; ++i) {
            if (AB[current][i] != 0 && !visited[i]) {
                next = i;
                break;
            }
        }

        // Si aucun sommet valide non visité n'est trouvé, choisir aléatoirement un sommet non visité et continuer la recherche des successeurs à partir de ce sommet
        if (next == -1) {
            for (int i = 0; i < length; ++i) {
                if (!visited[i]) {
                    next = i;
                    break;
                }
            }
        }

        // Si aucun sommet non visité n'est trouvé, revenir au sommet de départ
        if (next == -1) {
            next = 0; // Retourner au sommet de départ
        }

        current = next; // Définir le prochain sommet comme le sommet actuel

    } while (current != 0); // Continuer jusqu'à ce que nous revenions au sommet de départ

    // Suppression des sommets de valeur 0 dans le chromosome
    result.chromT.erase(std::remove(result.chromT.begin(), result.chromT.end(), 0), result.chromT.end());

    // Comparaison avec le parent 1 et ajout des clients manquants
    for (int i = 0; i < parent1.chromT.size(); ++i) {
        if (std::find(result.chromT.begin(), result.chromT.end(), parent1.chromT[i]) == result.chromT.end()) {
            result.chromT.push_back(parent1.chromT[i]); // Ajouter le client manquant au chromosome de l'enfant
        }
    }

    // Compléter l'individu avec l'algorithme Split
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverGrPX(Individual &result, const Individual &parent1, const Individual &parent2)
{
    int nbParents = 2;
    double routeCost[2][params.nbVehicles];
    double routeBenefit[2][params.nbVehicles];
    int tabAffectedClients[params.nbClients];

    for (int i=0; i<params.nbClients; i++){
        tabAffectedClients[i] = -1;
    }

    std::vector <Individual> vParents (2);

    vParents[0] = parent1;
    vParents[1] = parent2;


    for (int j=0; j<params.nbVehicles; j++){

        result.chromR[j].clear();
    }

    for (int k=0; k<nbParents; k++) {

        for (int l=0; l<params.nbVehicles; l++){
            routeCost[k][l]=0;
            routeBenefit[k][l]=0;
        }


        for (int l=0; l<params.nbVehicles; l++){


            int currentClient = 0;


            for (int customer : vParents[k].chromR[l]) {


                routeCost[k][l] += params.timeCost[currentClient][customer];
                routeBenefit[k][l] += params.timeCost[0][customer] * params.cli[customer].demand;



                currentClient = customer;




            }

            routeCost[k][l] += params.timeCost[0][currentClient];


        }


    }



    double valMax;
    int vehMax;

    int indice=0;

    for (int i=0; i<params.nbVehicles; i++) {



        Individual& currentParent = vParents[indice];
                double* currentRouteCost = routeCost[indice];
        double* currentRouteBenefit = routeBenefit[indice];

                valMax=-1;
                vehMax=(rand()/(double)RAND_MAX) * params.nbVehicles ;


        for (int j=0; j<params.nbVehicles; j++) {

            double currentVal;
            if(currentRouteCost[j] != 0){
                currentVal = (double)currentRouteBenefit[j]/ (double)currentRouteCost[j];
            } else {
                currentVal = -10;
            }


            if (currentVal>valMax) {
                valMax=currentVal;
                vehMax=j;
            }

        }



        for (int customer : currentParent.chromR[vehMax]) {


            if(tabAffectedClients[customer-1] < 0){

                tabAffectedClients[customer-1] = i;

                result.chromR[i].push_back(customer);

            }





        }


        routeCost[indice][vehMax]=0;
        routeBenefit[indice][vehMax]=0;


        int indice2;

        if(indice == 0){
            indice2 = 1;

        } else {
            indice2 = 0;
        }



        for (int l=0; l<params.nbVehicles; l++){
            routeCost[indice2][l]=0;
            routeBenefit[indice2][l]=0;
        }


        for (int l=0; l<params.nbVehicles; l++){


            int currentClient = 0;

            double totalDemande = 0;


            for (int customer : vParents[indice2].chromR[l]) {

                if(tabAffectedClients[customer-1] < 0){

                    routeCost[indice2][l] += params.timeCost[currentClient][customer] ;
                    routeBenefit[indice2][l] += params.timeCost[0][customer]* params.cli[customer].demand;

                    //routeBenefit[indice2][l] += 1;
                    //routeBenefit[indice2][l] += params.cli[customer].demand;

                    currentClient = customer;

                }


            }

            routeCost[indice2][l] += params.timeCost[0][currentClient];

        }

        if(indice == 0){
            indice = 1;

        } else {
            indice = 0;
        }
    }
    for (int i=0; i<params.nbClients; i++){

                if (tabAffectedClients[i]<0) {
                        int randVeh =(rand()/(double)RAND_MAX) * params.nbVehicles ;
            result.chromR[randVeh].push_back(i+1);
                }
        }

    result.evaluateCompleteCost(params);

}

std::vector<std::vector<int>> evalDistRoute(int r1, int r2, std::vector<int>  v1, std::vector<int>  v2){

    int c;

    int n = v1.size();
    int m = v2.size();


    std::vector<std::vector<int> > D (n+1, std::vector<int>(m+1));


    for (int i = 0; i < n+1; i++){
        D[i][0] = i;
    }


    for (int j = 0; j < m+1; j++){
        D[0][j] = j;
    }


    for (int i = 1; i < n+1; i++){

        for (int j = 1; j < m+1; j++){

            if (v1[i-1]  == v2[j-1]){

                c = 0;
            } else {

                c = 9999;

            }

            D[i][j] = std::min(std::min(D[i - 1][j] + 1, D[i][j - 1] + 1), D[i - 1][j - 1] + c);

        }

    }


    return D;
}


class Dmatrix {
    public:
        std::vector<std::vector<int> > tab;


};


void Genetic::crossoverPathRelinking(Individual & result, const Individual & parent1, const Individual & parent2)
{


        std::vector<std::vector<int> > A (params.nbVehicles, std::vector<int>(params.nbVehicles));

        std::vector <int> affectation = std::vector <int> (params.nbVehicles);
        std::vector <int> affectationBis = std::vector <int> (params.nbVehicles);

        std::vector<std::vector<int> > ifRoadReversed (params.nbVehicles, std::vector<int>(params.nbVehicles));



        std::vector<std::vector<Dmatrix>> storeDistMatrix(params.nbVehicles, std::vector<Dmatrix>(params.nbVehicles));




        for (int r1 = 0; r1 < params.nbVehicles; r1++){

            for (int r2 = 0; r2 < params.nbVehicles; r2++){


               std::vector<int> v1 = parent1.chromR[r1];
               std::vector<int> v2 = parent2.chromR[r2];

               int n = v1.size();
               int m = v2.size();


               std::vector<std::vector<int> > D1 = evalDistRoute(r1, r2, v1, v2);

               int d1 =  D1[n][m];

               std::reverse(v2.begin(),v2.end());



               std::vector<std::vector<int> > D2 = evalDistRoute(r1, r2, v1, v2);

               int d2 =  D2[n][m];



                if(d1 < d2){
                    A[r1][r2] = d1;
                    Dmatrix d;
                    d.tab = D1;

                    storeDistMatrix[r1][r2] = d;


                } else{
                    A[r1][r2] = d2;
                    ifRoadReversed[r1][r2] = 1;
                    Dmatrix d;
                    d.tab = D2;

                    storeDistMatrix[r1][r2] = d;
                }


            }

        }




        for (int c = 0; c < params.nbVehicles; c++){

            int minVal = 999999;
            int minI = -1;
            int minJ = -1;

            for (int i = 0; i < params.nbVehicles; i++){

                for (int j = 0; j < params.nbVehicles; j++){

                    if (A[i][j] < minVal){

                        minVal = A[i][j];
                        minI = i;
                        minJ = j;

                    }

                }

            }

            for (int k = 0; k < params.nbVehicles; k++){

                A[minI][k] = 9999;
                A[k][minJ] = 9999;

            }

            affectation[minI] = minJ;

            affectationBis[minJ] = minI;


        }


        for (int k = 0; k < params.nbVehicles; k++){
            int r1 = k;
            int r2 = int(affectation[k]);


        }


        std::vector<std::set <int>> toRemove (params.nbVehicles);

        std::vector<std::set <int>> toAdd (params.nbVehicles);

        std::vector<std::set <int>> F (params.nbVehicles);

        std::set <int> M;


        for (int k = 0; k < params.nbVehicles; k++){



            int cpt_move = 0;
            int r1 = k;
            int r2 = int(affectation[k]);

            std::vector<int> v1 = parent1.chromR[r1];
            std::vector<int> v2 = parent2.chromR[r2];

            int n = v1.size();
            int m = v2.size();



            if(ifRoadReversed[r1][r2] == 1){

//                 std::cout <<  "v2 Reversed " << std::endl ;

                std::reverse(v2.begin(),v2.end());

            }
            
            std::vector<std::vector<int> > D = storeDistMatrix[r1][r2].tab;


            int pos1 = n;
            int pos2 = m;


            int cpt = 0;




            while (pos1 != 0 || pos2!= 0){


                if(pos1 == 0){


                    toAdd[r2].insert(v2[pos2-1]);

                    M.insert(v2[pos2-1]);


                    pos2 = pos2 - 1;


                } else if(pos2 == 0) {

                    //M.push_back(v1[pos1-1]);

                    toRemove[r1].insert(v1[pos1-1]);

                    M.insert(v1[pos1-1]);

                   // M.push_back(0);

                    pos1 = pos1 - 1;


                } else if(D[pos1-1][pos2 - 1] == D[pos1][pos2]  && D[pos1-1][pos2 - 1] < D[pos1][pos2 - 1] && D[pos1-1][pos2-1] < D[pos1-1][pos2]){


                    F[r2].insert(v2[pos2-1]);

                    //F.push_back(0);
                    pos1 = pos1 - 1;
                    pos2 = pos2 - 1;




                } else if( D[pos1-1][pos2] < D[pos1-1][pos2 - 1]  && D[pos1-1][pos2] < D[pos1][pos2 - 1]){

                    //M.push_back(v1[pos1-1]);
                   // M.push_back(0);
                    toRemove[r1].insert(v1[pos1-1]);

                    M.insert(v1[pos1-1]);

                    pos1 = pos1 - 1;

                } else {

                    //M.push_back(v2[pos2-1]);

                    toAdd[r2].insert(v2[pos2-1]);

                    M.insert(v2[pos2-1]);


                    //M.push_back(0);
                    pos2 = pos2 - 1;

                }



            }



        }


        std::set< std::tuple<int,int,int>> relocateMoves;


        for (auto m : M) {
//             std::cout << "client " << m;

            for (int i = 0; i < params.nbVehicles; i++){

                if (toRemove[i].count(m)) {

                    for (int j = 0; j < params.nbVehicles; j++){

                        if (toAdd[j].count(m)) {

                            std::tuple <int, int, int> move = std::make_tuple(m, i, j);

                            relocateMoves.insert(move);
                        }

                    }

                }

            }

        }


        int nbRelocateMovesTotal = relocateMoves.size();
//         std::cout << "nbRelocateMovesTotal : " << nbRelocateMovesTotal << std::endl;

        int nbRelocateMovesTODO = nbRelocateMovesTotal/2;
//         std::cout << "nbRelocateMovesTODO : " << nbRelocateMovesTODO << std::endl;


        std::vector< std::tuple<int,int,int>> relocateMovesToPerform;

        std::sample(relocateMoves.begin(), relocateMoves.end(), std::back_inserter(relocateMovesToPerform), nbRelocateMovesTODO, std::mt19937{std::random_device{}()});


        std::random_shuffle(relocateMovesToPerform.begin(), relocateMovesToPerform.end());

        for (int k = 0; k < params.nbVehicles; k++){

            result.chromR[k] = parent1.chromR[k];

        }

        for (auto move : relocateMovesToPerform){

            int c = std::get<0>(move);
            int r1 = std::get<1>(move);
            int r2 = std::get<2>(move);



            result.chromR[r1].erase(std::remove(result.chromR[r1].begin(), result.chromR[r1].end(), c), result.chromR[r1].end());;



            std::vector<int> v2 = parent2.chromR[r2];

            if(ifRoadReversed[affectationBis[r2]][r2] == 1){

//                 std::cout <<  "v2 Reversed " << std::endl ;
                std::reverse(v2.begin(),v2.end());

            }


            int pivot = -1;

            bool notFound = true;

            int idx = 0;

            while(notFound){

                if (F[r2].count(v2[idx])) {

                    pivot = v2[idx];
                }

                if(v2[idx] == c){

                    notFound = false;
                }

                idx = idx+1;

            }

            F[r2].insert(c);

//             std::cout <<  "pivot " << pivot << std::endl ;


            notFound = true;
            idx = 0;

            if(pivot != -1){
                while(notFound){


                    if(result.chromR[affectationBis[r2]][idx] == pivot){

                        notFound = false;
                    }

                    idx = idx+1;

                }
            }

            result.chromR[affectationBis[r2]].insert(result.chromR[affectationBis[r2]].begin() + idx, c);



        }

        // Build up the rest of the Individual structure
        result.evaluateCompleteCost(params);

}

void Genetic::crossoverGrPX2(Individual &result, const Individual &parent1, const Individual &parent2)
{
    int nbParents = 2;
    double routeCost[2][params.nbVehicles];
    double routeBenefit[2][params.nbVehicles];
    int tabAffectedClients[params.nbClients];

    for (int i=0; i<params.nbClients; i++){
        tabAffectedClients[i] = -1;
    }

    std::vector <Individual> vParents (2);

    vParents[0] = parent1;
    vParents[1] = parent2;


    for (int j=0; j<params.nbVehicles; j++){

        result.chromR[j].clear();
    }

    for (int k=0; k<nbParents; k++) {

        for (int l=0; l<params.nbVehicles; l++){
            routeCost[k][l]=0;
            routeBenefit[k][l]=0;
        }


        for (int l=0; l<params.nbVehicles; l++){


            int currentClient = 0;


            for (int customer : vParents[k].chromR[l]) {


                routeCost[k][l] += params.timeCost[currentClient][customer];
                routeBenefit[k][l] += params.timeCost[0][customer] * params.cli[customer].demand;
                //routeBenefit[k][l] += params.timeCost[0][customer];

                currentClient = customer;




            }

            routeCost[k][l] += params.timeCost[0][currentClient];

        }


    }



    double valMax;
    int vehMax;

    int indice=0;

    int cpt = 0;

    for (int i=0; i<params.nbVehicles; i++) {



        Individual& currentParent = vParents[indice];
                double* currentRouteCost = routeCost[indice];
        double* currentRouteBenefit = routeBenefit[indice];

                valMax=-1;
                vehMax=(rand()/(double)RAND_MAX) * params.nbVehicles ;


        for (int j=0; j<params.nbVehicles; j++) {

            double currentVal;
            if(currentRouteCost[j] != 0){
                currentVal = (double)currentRouteBenefit[j]/ (double)currentRouteCost[j];
            } else {
                currentVal = -10;
            }


            if (currentVal>valMax) {
                valMax=currentVal;
                vehMax=j;
            }

        }



        for (int customer : currentParent.chromR[vehMax]) {


            if(tabAffectedClients[customer-1] < 0){

                tabAffectedClients[customer-1] = i;

                result.chromR[i].push_back(customer);

            }





        }


        routeCost[indice][vehMax]=0;
        routeBenefit[indice][vehMax]=0;


        int indice2;

        if(indice == 0){
            indice2 = 1;

        } else {
            indice2 = 0;
        }



        for (int l=0; l<params.nbVehicles; l++){
            routeCost[indice2][l]=0;
            routeBenefit[indice2][l]=0;
        }


        for (int l=0; l<params.nbVehicles; l++){


            int currentClient = 0;


            for (int customer : vParents[indice2].chromR[l]) {

                if(tabAffectedClients[customer-1] < 0){

                    routeCost[indice2][l] += params.timeCost[currentClient][customer];
                    routeBenefit[indice2][l] += params.timeCost[0][customer]* params.cli[customer].demand;
                    //routeBenefit[indice2][l] += params.timeCost[0][customer];
                    currentClient = customer;

                }


            }

            routeCost[indice2][l] += params.timeCost[0][currentClient];

        }

        if(indice == 0 and cpt == 1){
            indice = 1;

        } else if(indice == 0 and cpt == 0){
            cpt = 1;

        } else if(indice == 1){
            indice = 0;
            cpt = 0;
        }
    }


    for (int i=0; i<params.nbClients; i++){

                if (tabAffectedClients[i]<0) {
                        int randVeh =(rand()/(double)RAND_MAX) * params.nbVehicles ;
            result.chromR[randVeh].push_back(i+1);
                }
        }

    result.evaluateCompleteCost(params);

}

void Genetic::crossoverGOX(Individual & result, const Individual & parent1, const Individual & parent2) {
    std::vector<bool> freqClient(params.nbClients + 1, false);

    //posCut

    int pointCut = params.ap.nbCut;
    std::vector<int> cutPoints1(pointCut);  // Points de coupe pour parent1
    std::vector<int> cutPoints2(pointCut);  // Points de coupe pour parent2
    std::vector<std::vector<int>> segments1(pointCut);
    std::vector<std::vector<int>> segments2(pointCut);

    // Génération des points de coupe pour parent1
    if (params.ap.eqSeg == 1) {
        std::uniform_int_distribution<> distr(0, params.nbClients - 1);
        std::set<int> cuts1;
        while (cuts1.size() < pointCut) {
            cuts1.insert(distr(params.ran));
        }
        cutPoints1.assign(cuts1.begin(), cuts1.end());
    } else {
        int step = params.nbClients / pointCut;
        std::uniform_int_distribution<> distr(0, step - 1);
        int randomOffset = distr(params.ran);
        for (int i = 0; i < pointCut; ++i) {
            cutPoints1[i] = (randomOffset + i * step) % params.nbClients;
        }
    }
    std::sort(cutPoints1.begin(), cutPoints1.end());

    // Génération des points de coupe pour parent2
    if (params.ap.eqSeg == 1) {
        std::uniform_int_distribution<> distr(0, params.nbClients - 1);
        std::set<int> cuts2;
        while (cuts2.size() < pointCut) {
            cuts2.insert(distr(params.ran));
        }
        cutPoints2.assign(cuts2.begin(), cuts2.end());
    } else {
        int step = params.nbClients / pointCut;
        std::uniform_int_distribution<> distr(0, step - 1);
        int randomOffset = distr(params.ran);
        for (int i = 0; i < pointCut; ++i) {
            cutPoints2[i] = (randomOffset + i * step) % params.nbClients;
        }
    }
    std::sort(cutPoints2.begin(), cutPoints2.end());

    // Découpage de parent1
    for (int cutIndex = 0; cutIndex < cutPoints1.size(); ++cutIndex) {
        int start = cutPoints1[cutIndex];
        int end = (cutIndex + 1 < cutPoints1.size()) ? cutPoints1[cutIndex + 1] : (cutPoints1[0] + params.nbClients) % params.nbClients;

        for (int k = start; k != end; k = (k + 1) % params.nbClients) {
            segments1[cutIndex].push_back(parent1.chromT[k]);
        }
    }

    // Découpage de parent2
    for (int cutIndex = 0; cutIndex < cutPoints2.size(); ++cutIndex) {
        int start = cutPoints2[cutIndex];
        int end = (cutIndex + 1 < cutPoints2.size()) ? cutPoints2[cutIndex + 1] : (cutPoints2[0] + params.nbClients) % params.nbClients;

        for (int k = start; k != end; k = (k + 1) % params.nbClients) {
            segments2[cutIndex].push_back(parent2.chromT[k]);
        }
    }
    
    int useCostBenefit = params.ap.useCostBenefit;

    auto evaluateSegment = [&](const std::vector<int>& segment, std::vector<bool>& freqClient) {
        double totalCost = 0.0;
        double totalBenefit = 0.0;
        int currentClient = segment[0];
        int clientCount = 0; // Compteur de clients uniques

        for (int customer : segment) {
            if (!freqClient[customer]) {
                double cost = params.timeCost[currentClient][customer];
                double benefit = params.timeCost[0][customer] * params.cli[customer].demand;
                totalCost += cost;
                totalBenefit += benefit;
                currentClient = customer;
                clientCount++; // Incrémenter le compteur de clients uniques
            }
        }

        // Retourner le ratio Bénéfice/Coût si useCostBenefit est true, sinon retourner totalCost
        if (useCostBenefit == 0) {
            if (totalCost > 0) {
                return 1/totalCost; // Retourner le coût total (distance)
            }else{
                return std::numeric_limits<double>::min();
            }

        } else if (useCostBenefit == 1) {
            if (totalCost > 0) {
                return  totalBenefit/totalCost  ; // Ratio Bénéfice/Coût
            } else {
                return std::numeric_limits<double>::min(); // Coût zéro, retourner une valeur élevée
            }
        } else if (useCostBenefit == 2) {
            if (totalCost > 0) {
                return (totalBenefit - totalCost) / totalCost; // Ratio Bénéfice Net
            } else {
                return std::numeric_limits<double>::min(); // Coût zéro, retourner une valeur élevée
            }
        } else if (useCostBenefit == 3) {
            if (totalCost > 0) {
                return clientCount / totalCost; // Favorise les segments avec plus de clients et un coût minimal
            } else {
                return std::numeric_limits<double>::min();
            }
        }

        return 0.0; // Retourner 0 par défaut si aucune condition n'est satisfaite
    };

    int j = 0;
    int randomSegmentSelection = params.ap.randSelect;
    std::vector<std::vector<int>> chosenSegments;

    if (randomSegmentSelection == 0) {
        for (int i = 0; i < pointCut; ++i) {
            std::vector<int> currentSegment;

            if (i % 2 == 0) {
                for (int client : segments1[i]) {
                    if (!freqClient[client]) {
                        currentSegment.push_back(client);
                        freqClient[client] = true;
                    }
                }
            } else {
                for (int client : segments2[i]) {
                    if (!freqClient[client]) {
                        currentSegment.push_back(client);
                        freqClient[client] = true;
                    }
                }
            }

            if (!currentSegment.empty()) {
                chosenSegments.push_back(currentSegment);
            }
        }
    } else if (randomSegmentSelection == 1) {

        for (int i = 0; i < pointCut; ++i) {

            // Alternance entre les segments des deux parents
            std::vector<std::vector<int>>& segment = (i % 2 == 0) ? segments1 : segments2;

            double maxScore = std::numeric_limits<double>::min();
            int bestSegmentIndex = -1;

            // Parcours de tous les segments du parent sélectionné
            for (int segIndex = 0; segIndex < segment.size(); ++segIndex) {
                double score = evaluateSegment(segment[segIndex], freqClient);  // Évaluation du segment
                if (score > maxScore) {  // On cherche le segment avec le score maximum
                    maxScore = score;
                    bestSegmentIndex = segIndex;
                }
            }

            // Si un segment valide a été trouvé
            if (bestSegmentIndex != -1) {
                std::vector<int> currentSegment;
                // Ajout des clients non utilisés du meilleur segment
                for (int client : segment[bestSegmentIndex]) {
                    if (!freqClient[client]) {
                        currentSegment.push_back(client);
                        freqClient[client] = true;  // Marquer le client comme utilisé
                    }
                }

                // Si le segment sélectionné contient des clients, on l'ajoute
                if (!currentSegment.empty()) {
                    chosenSegments.push_back(currentSegment);
                }
            }
        }
    } else if (randomSegmentSelection == 2) {
        for (size_t i = 0; i < std::min(segments1.size(), segments2.size()); ++i) {
            double score1 = evaluateSegment(segments1[i], freqClient); // Évaluer le segment du parent 1
            double score2 = evaluateSegment(segments2[i], freqClient); // Évaluer le segment du parent 2

            bool chooseSegment1;

            // Si le critère est de minimiser le coût, un score plus bas est meilleur
            if (useCostBenefit == 0) {
                chooseSegment1 = score1 > score2; // Le score étant inversement proportionnel au coût, plus élevé est meilleur.
            } else {
                chooseSegment1 = score1 >= score2; // Pour les autres critères, un score plus élevé est toujours meilleur
            }

            std::vector<int> currentSegment;

            if (chooseSegment1) {
                // Ajouter le segment du parent 1
                for (int client : segments1[i]) {
                    if (!freqClient[client]) {
                        currentSegment.push_back(client);
                        freqClient[client] = true; // Marquer les clients comme utilisés
                    }
                }
            } else {
                // Ajouter le segment du parent 2
                for (int client : segments2[i]) {
                    if (!freqClient[client]) {
                        currentSegment.push_back(client);
                        freqClient[client] = true; // Marquer les clients comme utilisés
                    }
                }
            }

            if (!currentSegment.empty()) {
                chosenSegments.push_back(currentSegment); // Ajouter le segment choisi aux segments finaux
            }
        }
    }

    int insertionSegment = params.ap.insertSeg;

    if (insertionSegment == 0) {        
        result.chromT.clear();
        for (const auto& segment : chosenSegments) {
            for (int client : segment) {
                result.chromT.push_back(client);
            }
        }
    } else if (insertionSegment == 1) {
        std::shuffle(chosenSegments.begin(), chosenSegments.end(), std::mt19937{std::random_device{}()});
        result.chromT.clear();
        for (const auto& segment : chosenSegments) {
            for (int client : segment) {
                result.chromT.push_back(client);
            }
        }
    } else if (insertionSegment == 2) {

        // Initialisation
        result.chromT.clear();
        int currentClient = 0; // Début avec le client 0 (dépôt)

        // Boucle pour insérer les segments choisis de manière gloutonne
        while (!chosenSegments.empty()) {
            double minDist = std::numeric_limits<double>::max();
            int bestSegmentIndex = -1;

            for (size_t j = 0; j < chosenSegments.size(); ++j) {
                // Calculer la distance entre le dernier client ajouté et le premier client du segment
                double dist = params.timeCost[currentClient][chosenSegments[j][0]];

                // Vérifier si la distance est minimale
                if (dist < minDist) {
                    minDist = dist;
                    bestSegmentIndex = j;
                }
            }

            // Ajouter le meilleur segment trouvé
            if (bestSegmentIndex != -1) {
                const std::vector<int>& segment = chosenSegments[bestSegmentIndex];
                result.chromT.insert(result.chromT.end(), segment.begin(), segment.end());

                // Mettre à jour le client actuel
                currentClient = segment.back(); // Dernier client du segment ajouté

                // Supprimer le segment choisi de la liste
                chosenSegments.erase(chosenSegments.begin() + bestSegmentIndex);
            }
        }
    }

    for (int client : parent1.chromT) {
        if (!freqClient[client] && result.chromT.size() < params.nbClients) {
            result.chromT.push_back(client);
        }
    }

    split.generalSplit(result, parent1.eval.nbRoutes);
}

Genetic::Genetic(Params & params) :
        params(params),
        split(params),
        localSearch(params),
        population(params,this->split,this->localSearch),
        offspring(params){}
