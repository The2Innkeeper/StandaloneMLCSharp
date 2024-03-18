using System;
using System.Linq;
using System.Threading.Tasks;

class GeneticAlgorithm
{
    // Random number generator for various operations
    static Random rand = new Random();
    
    // Constants for the genetic algorithm
    const int PopulationSize = 100; // Size of the population
    const int Dimensions = 2; // Number of variables in the Ackley function (multivariable optimization)
    const int MaxGenerations = 100; // Maximum number of generations to run the algorithm
    const double MutationRate = 0.01; // Probability of mutation per gene
    const double CrossoverRate = 0.7; // Probability of crossover between pairs of individuals
    const double MinValue = -5.0; // Minimum value for each dimension in the Ackley function
    const double MaxValue = 5.0; // Maximum value for each dimension in the Ackley function

    static void Main(string[] args)
    {
        // Initialize population with random individuals
        double[][] population = new double[PopulationSize][];
        for (int i = 0; i < PopulationSize; i++)
        {
            population[i] = InitializeIndividual();
        }

        // Main loop of the genetic algorithm
        int generation = 0;
        while (generation < MaxGenerations)
        {
            // Evaluate fitness of each individual in parallel
            double[] fitness = new double[PopulationSize];
            Parallel.For(0, PopulationSize, i =>
            {
                fitness[i] = AckleyFunction(population[i]);
            });

            // Select individuals to create a mating pool
            double[][] selected = TournamentSelection(population, fitness);

            // Perform crossover to create new offspring
            double[][] offspring = Crossover(selected);

            // Apply mutation to the offspring
            Mutate(offspring);

            // Replacement with Elitism: keep the best individual from the current generation
            int bestIndex = Array.IndexOf(fitness, fitness.Min()); // Ackley function is minimized
            population[0] = population[bestIndex].ToArray(); // Copy the best individual
            for (int i = 1; i < PopulationSize; i++)
            {
                population[i] = offspring[i]; // Replace the rest of the population with offspring
            }

            // Print the best solution found in this generation
            Console.WriteLine($"Generation {generation}, Best Fitness: {fitness[bestIndex]:F2}, x: [{string.Join(", ", population[bestIndex].Select(x => x.ToString("F2")))}]");

            generation++;
        }
    }

    // Initialize an individual with random values within the specified range for each dimension
    static double[] InitializeIndividual()
    {
        return Enumerable.Range(0, Dimensions).Select(_ => rand.NextDouble() * (MaxValue - MinValue) + MinValue).ToArray();
    }

    // The Ackley function, a complex mathematical function used as a benchmark in optimization
    static double AckleyFunction(double[] x)
    {
        double a = 20, b = 0.2, c = 2 * Math.PI;
        double sum1 = 0.0, sum2 = 0.0;
        foreach (var xi in x)
        {
            sum1 += xi * xi; // Sum of squares
            sum2 += Math.Cos(c * xi); // Sum of cosines
        }
        double d = x.Length;
        // The Ackley function formula
        return -a * Math.Exp(-b * Math.Sqrt(sum1 / d)) - Math.Exp(sum2 / d) + a + Math.E;
    }

    // Tournament selection method for selecting individuals based on their fitness
    static double[][] TournamentSelection(double[][] population, double[] fitness)
    {
        double[][] selected = new double[PopulationSize][];
        for (int i = 0; i < PopulationSize; i++)
        {
            // Randomly select two individuals and choose the better one
            int first = rand.Next(PopulationSize);
            int second = rand.Next(PopulationSize);
            selected[i] = fitness[first] < fitness[second] ? population[first] : population[second];
        }
        return selected;
    }

    // Crossover method to combine genetic information from pairs of individuals
    static double[][] Crossover(double[][] selected)
    {
        double[][] offspring = new double[PopulationSize][];
        for (int i = 0; i < PopulationSize; i += 2)
        {
            offspring[i] = new double[Dimensions];
            offspring[i + 1] = new double[Dimensions];
            for (int j = 0; j < Dimensions; j++)
            {
                // Perform crossover with a certain probability
                if (rand.NextDouble() < CrossoverRate)
                {
                    offspring[i][j] = selected[i][j];
                    offspring[i + 1][j] = selected[i + 1][j];
                }
                else
                {
                    // Swap genes if crossover does not occur
                    offspring[i][j] = selected[i + 1][j];
                    offspring[i + 1][j] = selected[i][j];
                }
            }
        }
        return offspring;
    }

    // Mutation method to introduce random changes to the offspring
    static void Mutate(double[][] offspring)
    {
        for (int i = 0; i < offspring.Length; i++)
        {
            for (int j = 0; j < Dimensions; j++)
            {
                // With a certain probability, mutate a gene
                if (rand.NextDouble() < MutationRate)
                {
                    // Apply a small change to the gene value
                    offspring[i][j] += (rand.NextDouble() * 2 - 1) * (MaxValue - MinValue) * MutationRate;
                    // Ensure the mutated gene is within the predefined bounds
                    offspring[i][j] = Math.Max(MinValue, Math.Min(MaxValue, offspring[i][j]));
                }
            }
        }
    }
}
