package bio

val bases: Array<Char> = arrayOf('A', 'C', 'G', 'T')
val generation: Array<String> = arrayOf("AAAA", "AAAa", "AAaa", "AaAa", "Aaaa", "aaaa")

/**
 * FASTA formatted collection of DNA.
 */
open class FASTA<D, T>(list: Map<D, T>) {
    val list = list
}

open class Chemical<T>(strain:Array<T>) {
    val strain = strain
}
class Protein(strain: Array<Char>): Chemical<Char>(strain) {
    val strainString = strain.toString()

    fun motif(): List<Int> {
        var index = 0
        val r = Regex("N[^P][S|T][^P]")
        var res = listOf<Int>()
        while (r.find(strainString, index) != null) {
            res += r.find(strainString, index)!!.range.first + 1
            index = r.find(strainString, index)!!.range.first + 1
        }
        return res
    }
}

class Protein_FASTA(list: Map<String, Protein>) : FASTA<String, Protein>(list) {
}

class DNA_FASTA(list: Map<String, DNA>) : FASTA<String, DNA>(list) {
    /**
     * Returns profile of DNA strains given in FASTA format
     */
    fun profile(): Map<Char, Array<Int>> {
        val len: Int = list.values.first().strain.length
        var profile: Map<Char, Array<Int>> = hashMapOf('A' to Array(len, { 0 }), 'C' to Array(len, { 0 }), 'G' to Array(len, { 0 }), 'T' to Array(len, { 0 }))
        for (i in 0..len - 1)
            for (item in list.values) {
                var tmp: Char = item.strain[i]
                profile[tmp]!![i] += 1
            }
        return profile
    }

    /**
     * Returns consensus of DNA strains given in FASTA format
     */
    fun consensus(): String {
        val profile = profile()
        val len: Int = list.values.first().strain.length
        var res: String = ""
        for (i in 0..len - 1) {
            var max: Int = 0
            var maxs = 'A'
            for (base in bases) {
                if (profile[base]!![i] > max) {
                    max = profile[base]!![i]
                    maxs = base
                }
            }
            res += maxs
        }
        return (res)
    }

    /**
     * Returns map of all possible prefix/suffix pairs in given DNA strains set
     */
    fun overlap(): Map<String, String> {
        var res: Map<String, String> = hashMapOf()
        var pre: Map<String, DNA>
        var suf: Map<String, DNA>
        for (pre in list.keys)
            for (suf in list.keys)
                if ((list[pre]!!.strain.takeLast(3) == list[suf]!!.strain.take(3)) and (pre != suf))
                    res += Pair(pre, suf)
        return res
    }

    /**
     * Returns one of the motives shared by all strains in the DNA strains set
     */
    fun shared_motif(): String {
        var data = list.values.toTypedArray()
        var res: String = ""
        if ((data.size > 1) and (data[0].strain.length > 0))
            for (i in 0..data[0].strain.length - 1)
                for (j in 0..data[0].strain.length - i) {
                    if ((j > res.length) and (data.all { it -> it.strain.contains(data[0].strain.drop(i).take(j)) })) {
                        res = data[0].strain.drop(i).take(j)
                    }
                }
        return res
    }
}

/**
 * DNA in a form of the string of bases.
 */
class DNA(strain: String) {
    val complements: Map<Char, Char> = hashMapOf('A' to 'T', 'C' to 'G', 'G' to 'C', 'T' to 'A')
    val strain: String = strain
    val transcription: RNA = RNA(strain.replace('T', 'U'))
    val gccontent: Float = ((nucleotidesCount()['C'] as Int + nucleotidesCount()['G'] as Int).toFloat() / strain.length.toFloat()) * 100

    fun nucleotidesCount(): Map<Char, Int> {
        //TODO: resolve the issue with nullability and typecasting
        var res: Map<Char, Int> = mapOf()
        for (item in strain.groupBy { it }.toSortedMap()) {
            res += Pair(item.key, item.value.size)
        }
        for (key in bases) {
            if (res[key] == null)
                res += Pair(key, 0)
        }
        return res
    }

    /**
     * Returns complementary DNA
     */
    fun complement(): DNA {
        var res: String = ""
        for (base in strain)
            res += complements[base]
        return DNA(res.reversed())
    }

    /**
     * Returns hammingDistance between this DNA and given
     */
    fun hammingDistance(t: DNA): Int {
        var res: Int = 0
        for ((a, b) in strain.zip(t.strain)) {
            if (a != b)
                res++
        }
        return res
    }

    /**
     * Returns motif t occurences in this DNA
     */
    fun motif(t: DNA): Array<Int> {
        var res: Array<Int> = arrayOf()
        for (i in 0..strain.length - 1 - t.strain.length)
            if (strain.substring(i, i + t.strain.length) == t.strain)
                res += i + 1
        return res
    }
}

/**
 * RNA in a form of the string of bases
 */
class RNA(strain: String) {
    val strain: String = strain
    val codons: Map<String, Char> = hashMapOf("UUU" to 'F', "CUU" to 'L', "AUU" to 'I', "GUU" to 'V', "UUC" to 'F', "CUC" to 'L', "AUC" to 'I', "GUC" to 'V', "UUA" to 'L',
            "CUA" to 'L', "AUA" to 'I', "GUA" to 'V', "UUG" to 'L', "CUG" to 'L', "AUG" to 'M', "GUG" to 'V', "UCU" to 'S', "CCU" to 'P',
            "ACU" to 'T', "GCU" to 'A', "UCC" to 'S', "CCC" to 'P', "ACC" to 'T', "GCC" to 'A', "UCA" to 'S', "CCA" to 'P', "ACA" to 'T',
            "GCA" to 'A', "UCG" to 'S', "CCG" to 'P', "ACG" to 'T', "GCG" to 'A', "UAU" to 'Y', "CAU" to 'H', "AAU" to 'N', "GAU" to 'D',
            "UAC" to 'Y', "CAC" to 'H', "AAC" to 'N', "GAC" to 'D', "UAA" to '!', "CAA" to 'Q', "AAA" to 'K', "GAA" to 'E', "UAG" to '!',
            "CAG" to 'Q', "AAG" to 'K', "GAG" to 'E', "UGU" to 'C', "CGU" to 'R', "AGU" to 'S', "GGU" to 'G', "UGC" to 'C', "CGC" to 'R',
            "AGC" to 'S', "GGC" to 'G', "UGA" to '!', "CGA" to 'R', "AGA" to 'R', "GGA" to 'G', "UGG" to 'W', "CGG" to 'R', "AGG" to 'R',
            "GGG" to 'G')

    fun protein(): String {
        var res: String = ""
        for (i in 0..strain.length / 3 - 1)
            res += codons[strain[i * 3] + strain[i * 3 + 1].toString() + strain[i * 3 + 2]]
        return res
    }
}

/**
 * Return number of rabbits after n months if each pair spawns k more after 2 months of living
 */
fun rabbits(n: Int, k: Int): Long {
    when (n) {
        1, 2 -> return 1L
        else -> return rabbits(n - 1, k) + rabbits(n - 2, k) * k
    }
}

/**
 * Return number of rabbits after n months if each pair spawns 1 more after 2 months of living and die after m months
 */
fun mortalrabbits(n: Int, m: Int = 1): Long {
    var ages: Array<Long> = Array(m, { 0L })
    ages[0] = 1
    for (i in 1..n - 1) {
        var tmp: Long = ages.drop(1).sum()
        for (j in (m - 1).downTo(1))
            ages[j] = ages[j - 1]
        ages[0] = tmp
    }
    return ages.sum()
}

/**
 * Return possibility of the birth of the specimen carrying dominant allele
 */
fun mendel(k: Int, m: Int, n: Int): Double {
    var k = k.toDouble()
    var m = k.toDouble()
    var n = n.toDouble()
    return (k * (k - 1 + 2 * m + 2 * n) + m * (0.75 * (m - 1) + n)) * 1.0 / ((k + m + n) * (k + m + n - 1) * 1.0)
}

/**
 * Returns expected number of offspring displaying the dominant phenotype in the next generation, under the assumption that every couple has exactly two offspring.
 */
fun offspring(AAAA: Int, AAAa: Int, AAaa: Int, AaAa: Int, Aaaa: Int, aaaa: Int): Double {
    return (2.0 * ((AAAA + AAAa + AAaa).toDouble() + AaAa.toDouble() * 0.75 + Aaaa.toDouble() * 0.5))
}

/**
 * Tail-recursive factorial
 */
tailrec fun factorial(n: Int, offset: Int = 1): Int {
    if (n == 0)
        return offset
    else
        return factorial(n - 1, (offset * n))
}

fun power(n: Int, p: Int): Int {
    var res: Int = 1
    for (i in 1..p) {
        res = res * n
    }
    return res
}

fun power(n: Double, p: Int): Double {
    var res: Double = 1.0
    for (i in 1..p) {
        res = res * n
    }
    return res
}

/**
 * TODO: Write comment about combinatorics sense of these
 */
fun binomial(n: Int, k: Int): Int {
    return factorial(n) / (factorial(k) * factorial(n - k))
}


/**
 * Bernulli probability
 */
fun probability(n: Int, k: Int, p: Double): Double {
    return binomial(power(2, k), n) * power(p, n) * power((1 - p), (power(2, k) - n))
}

/**
 * TODO:independentAlleles, works incorrectly now
 */
fun independentAlleles(n: Int, k: Int): Double {
    var sum: Double = 0.0
    for (i in 0..n - 1)
        sum += probability(i, k, 0.25)
    return 1 - sum
}