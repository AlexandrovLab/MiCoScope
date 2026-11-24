/**
 * Which MetaPhlAn database, index and environment this launch uses.
 *
 * The default is the newest database the site configures (--metaphlan_db; vJan26 on TSCC),
 * run with MetaPhlAn 4.2. HUMAnN 3.9 refuses profiles from databases newer than vJun23
 * ("not generated with the database version v3 or vJun23"; probed 23 Nov 2026), so a launch
 * that asks for functional profiles (--skip_humann3 false) with --metaphlan_index auto
 * switches to --metaphlan_db_vjun23 and MetaPhlAn 4.1.1. Explicit settings always win:
 * an explicit --metaphlan_index is used as given (HUMAnN is then skipped with a warning when
 * it is incompatible; main.nf humannDbGuard), and an explicit --metaphlan4_env is used as given.
 *
 * --metaphlan_index auto takes the index name from the single mpa_*.pkl in the chosen
 * database directory, so a run pointed at an older database keeps working unchanged.
 * Every method is pure over params and the database directory listing; modules and main.nf
 * call the same methods, so the launch check, the tasks and the run record cannot disagree.
 */
class Mpa {

    /** --skip_humann3 false, in either the Boolean (24.10) or String (26.x) form. */
    static boolean humannRequested(Map params) {
        def v = params.skip_humann3
        return v == false || v?.toString() == 'false'
    }

    /** Profiles from this index are accepted by HUMAnN 3.9. */
    static boolean humannCompatible(String index) {
        return index != null && (index.contains('vJun23') || index.startsWith('mpa_v3'))
    }

    static boolean autoIndex(Map params) {
        return params.metaphlan_index == null || params.metaphlan_index.toString() == 'auto'
    }

    /** The single mpa_*.pkl in dir, without .pkl; null when there is none or more than one. */
    static String detect(Object dir) {
        if (dir == null) return null
        def d = new File(dir.toString())
        if (!d.isDirectory()) return null
        def pkls = (d.listFiles()?.toList() ?: []).findAll { it.isFile() && it.name ==~ /mpa_.*\.pkl/ }
        return pkls.size() == 1 ? pkls[0].name[0..-5] : null
    }

    /** True when this launch moves to the vJun23 database so HUMAnN can read the profiles. */
    static boolean switched(Map params) {
        return autoIndex(params) && humannRequested(params) &&
               !humannCompatible(detect(params.metaphlan_db))
    }

    /** The database directory MetaPhlAn (and HUMAnN's MetaPhlAn options) read. */
    static String dbDir(Map params) {
        def dir = switched(params) ? params.metaphlan_db_vjun23 : params.metaphlan_db
        return dir?.toString()
    }

    /** The index name inside dbDir(). */
    static String index(Map params) {
        return autoIndex(params) ? detect(dbDir(params)) : params.metaphlan_index.toString()
    }

    /** The conda env spec: 4.1.1 for vJun23-era indexes, 4.2 otherwise, unless set explicitly. */
    static String env(Map params) {
        def explicit = params.metaphlan4_env?.toString()
        if (explicit && explicit != 'auto') return explicit
        return humannCompatible(index(params)) ? params.metaphlan41_env.toString()
                                               : params.metaphlan42_env.toString()
    }

    /** null when consistent, else the message main.nf stops the launch with. */
    static String check(Map params) {
        if (switched(params) && !params.metaphlan_db_vjun23)
            return "--skip_humann3 false needs a vJun23-era MetaPhlAn database for HUMAnN 3.9: set " +
                   "--metaphlan_db_vjun23 (the database in --metaphlan_db, '${params.metaphlan_db}', " +
                   "is newer), or pass --metaphlan_index to choose explicitly"
        def dir = dbDir(params)
        def idx = index(params)
        if (idx == null)
            return "--metaphlan_index auto: ${dir} must hold exactly one mpa_*.pkl; " +
                   "pass --metaphlan_index to choose one"
        if (new File(dir).isDirectory() && !new File(dir, "${idx}.pkl").exists())
            return "MetaPhlAn database ${dir} has no ${idx}.pkl (check --metaphlan_index)"
        return null
    }

    /** One line for the launch log. */
    static String describe(Map params) {
        def why = switched(params) ? ' (vJun23: HUMAnN requested)' : ''
        return "MetaPhlAn: ${index(params)} in ${dbDir(params)}${why}"
    }
}
