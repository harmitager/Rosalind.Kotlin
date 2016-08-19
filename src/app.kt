import bio.*
import java.io.File
import java.net.URL

fun main(args: Array<String>) {
    var input = File("C:\\Users\\araspopov\\Dropbox\\Rosalind\\src\\input.txt").readLines()
    var input2 = listOf<String>()
    for (url in input)
        for (item in URL("http://www.uniprot.org/uniprot/" + url + ".fasta").readText().split('\n'))
            if (item != "")
                input2 += item
    var key: String = ""
    var value: String = ""
    var L: Map<String, DNA> = mapOf()
    for (item in input2) {
        if (item[0] == '>') {
            if (value.length != 0) {
                L = L.plus(Pair(key, DNA(value)))
                value = ""
            }
            key = item
        } else {
            value += item
        }
    }
    L = L.plus(Pair(key, DNA(value)))
}