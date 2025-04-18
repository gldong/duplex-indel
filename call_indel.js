#!/usr/bin/env k8

/************
 * Performs double-stranded and single-stranded Indel calling from Tn5 duplex sequencing data. 
 * 
 * NOTE: 
 * This code is inspired by and adapted from https://github.com/lh3/lianti (21a15c8) which was
 * originally written by Heng Li to perform duplex SNV calling. We thank Heng Li for allowing us 
 * to adapt this code for indel calling in this repo.
 ************/

var version = "r158";

/************
 * getopt() *
 ************/

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

/*************************************
 * Parameters & command-line parsing *
 *************************************/

var c, min_mapq = 50, flt_win = 100, n_bulk = 1, is_hap_cell = false, show_flt = false, auto_only = false;
var min_dp_alt_cell = 5, min_dp_alt_strand_cell = 2, min_ab_cell = 0.2, max_lt_cell = 1, min_end_len = 10, min_joint_cell = 2, min_joint_strand_cell = 1;
var min_dp_bulk = 20, min_het_dp_bulk = 8, max_alt_dp_bulk = 0, min_het_ab_bulk = 0.3, is_hap_bulk = false;
var min_dp_dmg_strand = 4;
var fn_var = null, fn_hap = null, fn_excl = null, fn_rep = null;
var max_dp_ref_cell = 0; //number of REF reads allowed at NV site

while ((c = getopt(arguments, "h:A:b:v:D:e:Hl:a:s:w:m:Fr:uL:B:S:Pj:J:R:T:")) != null) {
	if (c == 'b') n_bulk = parseInt(getopt.arg);
	else if (c == 'H') is_hap_cell = true;
	else if (c == 'h') fn_hap = getopt.arg;
	else if (c == 'e') fn_excl = getopt.arg;
	else if (c == 'v') fn_var = getopt.arg;
	else if (c == 'r') fn_rep = getopt.arg;
	else if (c == 'F') show_flt = true;
	else if (c == 'u') auto_only = true;
	else if (c == 'a') min_dp_alt_cell = parseInt(getopt.arg);
	else if (c == 's') min_dp_alt_strand_cell = parseInt(getopt.arg);
	else if (c == 'w') flt_win = parseInt(getopt.arg);
	else if (c == 'S') min_dp_dmg_strand = parseInt(getopt.arg);
	else if (c == 'B') min_ab_cell = parseFloat(getopt.arg);
	else if (c == 'l') max_lt_cell = parseInt(getopt.arg);
	else if (c == 'L') min_end_len = parseInt(getopt.arg);
	else if (c == 'D') min_dp_bulk = parseInt(getopt.arg);
	else if (c == 'A') min_het_dp_bulk = parseInt(getopt.arg);
	else if (c == 'm') max_alt_dp_bulk = parseInt(getopt.arg);
	else if (c == 'P') is_hap_bulk = is_hap_cell = true;
	else if (c == 'j') min_joint_cell = parseInt(getopt.arg);
	else if (c == 'J') min_joint_strand_cell = parseInt(getopt.arg);
	else if (c == 'R') max_dp_ref_cell = parseInt(getopt.arg);
	else if (c == 'T') flt_lv = parseInt(getopt.arg);
}

if (min_dp_alt_strand_cell * 2 > min_dp_alt_cell)
	throw("2 * {-s} should not be larger than {-a}");

if (arguments.length - getopt.ind == 0) {
	print("Usage: plp-joint.js [options] <joint.vcf>");
	print("Options:");
	print("  General:");
	print("    -b INT    number of bulk samples [1]");
	print("    -h FILE   samples in FILE are haploid []");
	print("    -H        mark all single-cell samples as haploid");
	print("    -e FILE   exclude samples contained in FILE []");
	print("    -v FILE   exclude positions in VCF FILE []");
	print("    -r FILE   cell replicates []");
	print("    -F        print Indels filtered by -w and -v");
	print("    -u        process autosomes only");
	print("  Cell:");
	print("    -a INT    min ALT read depth to call an Indel [" + min_dp_alt_cell + "]");
	print("    -s INT    min ALT read depth per strand [" + min_dp_alt_strand_cell + "]");
	print("    -l INT    max LIANTI conflicting reads [" + max_lt_cell + "]");
	print("    -L INT    min distance towards the end of a read [" + min_end_len + "]");
	print("    -w INT    size of window to filter clustered Indels [" + flt_win + "]");
	print("    -S INT    min strand depth at candidate DNA damages [" + min_dp_dmg_strand + "]");
	print("    -B FLOAT  min ALT allele balance [" + min_ab_cell + "]");
	print("    -j INT    min allele depth to call joint Indels [" + min_joint_cell + "]");
	print("    -J INT    min allele depth on both strands to call joint Indels [" + min_joint_strand_cell + "]");
	print("    -R INT    max REF read depth to call an Indel [" + max_dp_ref_cell + "]");
	print("    -T INT    filtering stringency level from 1 (most stringent) to 3 (least stringent) [" + flt_lv + "]");
	print("  Bulk:");
	print("    -D INT    min bulk read depth [" + min_dp_bulk + "]");
	print("    -A INT    min bulk ALT read depth to call a het [" + min_het_dp_bulk + "]");
	print("    -m INT    max bulk ALT read depth to call an Indel [" + max_alt_dp_bulk + "]");
	print("    -P        the bulk is haploid");
	exit(1);
}

print('CL', 'plp-joint.js ' + arguments.join(" "));
print('VN', version);
print('CC', 'SM  sample name (each sample)');
print('CC', 'NV  somatic Indels');
print('CC', 'NN  number of called somatic Indels (each)');
print('CC', 'NR  false negative rate (each)');
print('CC', 'NC  number of somatic Indels after FNR correction (each)');
print('CC', 'NA  alignment somatic Indels');
print('CC', 'DV  DNA damages or amplification errors');
print('CC', 'DN  number of called damages/errors (each)');
print('CC', 'DR  false negative rate of damages/errors (each)');
print('CC', 'DC  number of damages/errors after FNR correction (each)');
print('CC');

/***********************
 * Auxiliary functions *
 ***********************/

function read_list(fn)
{
	if (fn == null || fn == "") return {};
	var buf = new Bytes();
	var file = fn == '-'? new File() : new File(fn);
	var h = {};
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		h[t[0]] = 1;
	}
	file.close();
	buf.destroy();
	return h;
}

function aggregate_calls(x, cell_meta, is_hap_bulk)
{
	var bulk_ad = [0, 0], bulk_alt = [], cell_hit_jv = [], cell_hit_nv = [];
	for (var i = 0; i < x.bulk.length; ++i)
		bulk_ad[0] += x.bulk[i].ad[0], bulk_ad[1] += x.bulk[i].ad[1], bulk_alt.push(x.bulk[i].ad[1]);
	if (bulk_ad[1] != 0) bulk_ad[1] = bulk_alt.join(":");
	for (var i = 0; i < x.cell.length; ++i) {
		var c = x.cell[i];
		if (!x.flt && x.n_joint_alt >= 2) {
			var b;
			if (c.flt || c.dp == 0) b = '.';
			else if (c.ad[0] > 0 && c.ad[1] >= min_joint_cell) b = '1';
			else if (c.ad[1] >= min_joint_cell) b = is_hap_bulk? '6' : cell_meta[i].ploidy == 1? '4' : '2';
			else if (c.ad[1] > 0) b = '.';
			else b = is_hap_bulk? '5' : cell_meta[i].ploidy == 1? '3' : '0';
			cell_meta[i].calls.push(b);
			if (c.ad[1] >= min_joint_cell)
				cell_hit_jv.push([cell_meta[i].name, c.adf[1], c.adr[1], c.count_ltpos, c.uniq_ltpos, c.min_mg_start, c.max_mg_end].join(":")); //Tn5 sites & merged window
		}
		if (!c.flt && c.alt && !x.flt) {
			++cell_meta[i].indel;
			cell_hit_nv.push([cell_meta[i].name, c.adf[1], c.adr[1], c.count_ltpos, c.uniq_ltpos, c.min_mg_start, c.max_mg_end].join(":")); //Tn5 sites & merged window
		}
	}
	if (cell_hit_nv.length > 0)
		print('NV', x.ctg, x.pos, x.ref, x.alt, bulk_ad.join("\t"), cell_hit_nv.length, cell_hit_nv.join(","));
	if (cell_hit_jv.length > 0)
		print('JV', x.ctg, x.pos, x.ref, x.alt, bulk_ad.join("\t"), cell_hit_jv.length, cell_hit_jv.join(","));
}

/********
 * Main *
 ********/

var file, buf = new Bytes(), re_auto = new RegExp('^(chr)?([0-9]+)$');

var var_map = new Map();
if (fn_var != null) {
	warn('Reading sites to filter...');
	file = new File(fn_var);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t[0][0] == '#') continue;
		var_map.put(t[0] + ':' + t[1]);
	}
	file.close();
}

var rep_str = {};
if (fn_rep != null) {
	file = new File(fn_rep);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split(/\s+/);
		for (var i = 1; i < t.length; ++i)
			rep_str[t[i]] = t[0];
	}
	file.close();
}

var sample_excl = read_list(fn_excl);
var sample_hap = read_list(fn_hap);
var col2cell = [];
var cell_meta = [];

warn('Calling...');
file = arguments[getopt.ind] == '-'? new File() : new File(arguments[getopt.ind]);
var rep_id = [], last = [], last_bulk = [], n_het_bulk = 0, n_hom_bulk = 0, n_het_bulk_detected = 0;
while (file.readline(buf) >= 0) {
	var m, t = buf.toString().split("\t");
	if (t[0] == '#CHROM') { // parse the sample line
		var sample_name = [];
		for (var i = 9 + n_bulk; i < t.length; ++i) {
			var s1 = t[i], s2 = s1.replace(/\.bam$/, "");
			if (sample_excl[s1] || sample_excl[s2]) continue;
			if (rep_str[s1] || rep_str[s2]) continue;
			var pl = is_hap_cell || sample_hap[s1] || sample_hap[s2]? 1 : 2;
			cell_meta.push({ name:s2, ploidy:pl, col:i, ado:[0,0], fn:0, indel:0, dmg:0, dmg_fp:0, dmg_fn:[0, 0], calls:[] });
			sample_name.push(s2); // for printing only
		}
		for (var i = 0; i < cell_meta.length; ++i)
			col2cell[cell_meta[i].col] = i;
		// construct rep_id
		for (var i = 0; i < t.length; ++i)
			rep_id[i] = i;
		if (fn_rep != null) {
			var sample2id = {};
			for (var i = 9; i < t.length; ++i) {
				var s1 = t[i], s2 = s1.replace(/\.bam$/, "");
				sample2id[s1] = sample2id[s2] = i;
			}
			for (var i = 9; i < t.length; ++i) {
				var s1 = t[i], s2 = s1.replace(/\.bam$/, "");
				if (rep_str[s1] || rep_str[s2]) {
					var s3 = rep_str[s1] != null? rep_str[s1] : rep_str[s2];
					if (sample2id[s3] != null) rep_id[i] = sample2id[s3];
				}
			}
		}
		print('SM', sample_name.join("\t"));
		continue;
	} else if (t[0][0] == '#') continue; // skip header

	if (auto_only && !re_auto.test(t[0])) continue;
	t[1] = parseInt(t[1]);

	var flt_bulks = false, flt_indel = false;

	// skip bad sites: mapQ
	if ((m = /AMQ=([\d,]+)/.exec(t[7])) != null) {
		var s = m[1].split(","), flt = false;
		for (var j = 0; j < s.length; ++j)
			if (parseInt(s[j]) < min_mapq)
				flt = true;
		if (flt) flt_bulks = flt_indel = true;
	}

	// parse the FORMAT field
	var fmt = t[8].split(":"), fmt_hash = {};
	for (var i = 0; i < fmt.length; ++i)
		fmt_hash[fmt[i]] = i;
	var fmt_ltdrop = fmt_hash["LTDROP"];
	var fmt_alen = fmt_hash["ALEN"];
	var fmt_adf = fmt_hash["ADF"];
	var fmt_adr = fmt_hash["ADR"];
	if (fmt_adf == null || fmt_adr == null)
		throw Error('missing ADF or ADR in FORMAT');
	var fmt_ltpos = fmt_hash["LTPOS"]; //Tn5 sites
	var fmt_mgpos = fmt_hash["MGPOS"]; //merged window

	// parse VCF (this part works with multiple ALT alleles)
	var cell = [], bulk = [];
	for (var i = 9; i < t.length; ++i) {
		var cell_id = col2cell[rep_id[i]];
 		if (i >= 9 + n_bulk && cell_id == null) continue; // exclude this sample
		var s = t[i].split(":");
		var lt = fmt_ltdrop != null && s[fmt_ltdrop] != '.'? parseInt(s[fmt_ltdrop]) : 0;
		var adf = s[fmt_adf].split(",");
		var adr = s[fmt_adr].split(",");
		var ad = [], dp = 0;
		if (adf.length != adr.length) throw Error("Inconsistent VCF");
		var dp_ref = 0, dp_alt = 0;
		for (var j = 0; j < adf.length; ++j) {
			adf[j] = parseInt(adf[j]);
			adr[j] = parseInt(adr[j]);
			if (j == 0) dp_ref += adf[j] + adr[j];
			else dp_alt += adf[j] + adr[j];
			ad[j] = adf[j] + adr[j];
			dp += ad[j];
		}
		if (i < 9 + n_bulk) {
			bulk.push({ dp:dp, ad:ad, adf:adf, adr:adr });
		} else {
			var flt = false, flt_dmg = false;
			if (cell_meta[cell_id].ploidy == 1 && dp_alt > 0 && dp_ref > 0) flt = true; // two alleles in a haploid cell; flt_dmg is not affected by this
			if (lt > max_lt_cell) flt = true;
			if (fmt_alen != null && s[fmt_alen] != '.') {
				var u = s[fmt_alen].split(",");
				for (var j = 1; j < u.length; ++j)
					if (u[j] != '.' && parseFloat(u[j]) < min_end_len)
						flt = flt_dmg = true;
			}
			//left and right ends of read/read-pair
			var ltpos = [], obj_count = {}, uniq_ltpos = [], count_ltpos = 0;
			if (fmt_ltpos != null && s[fmt_ltpos] != '.' && s[fmt_ltpos] != ""){
				ltpos = s[fmt_ltpos].split(",");
				for (var j = 0; j < ltpos.length; ++j) {
					var key = ltpos[j];
					if (key == "" || key == "." || key == "|") continue;
					if (obj_count.hasOwnProperty(key)) obj_count[key]++;
					else obj_count[key] = 1;
				}
				uniq_ltpos = Object.getOwnPropertyNames(obj_count);
				count_ltpos = uniq_ltpos.length;
				uniq_ltpos = uniq_ltpos.toString();
			}
			//left and right sites of merged window
			var mgpos = [], min_mg_start = Infinity, max_mg_end = 0;
			if (fmt_mgpos != null && s[fmt_mgpos] != '.' && s[fmt_mgpos] != ""){
				mgpos = s[fmt_mgpos].split(",");
				for (var j = 0; j < mgpos.length; ++j) {
					var key = mgpos[j];
					if (key == "" || key == "." || key == "|") continue;
					var mg_start = key.split("|")[0]; 
					var mg_end = key.split("|")[1];
					mg_start = parseInt(mg_start);
					mg_end = parseInt(mg_end);
					if (mg_start < min_mg_start) min_mg_start = mg_start;
					if (mg_end > max_mg_end) max_mg_end = mg_end;
				}
				var is_deletion = true, indel_start = t[1] + 1, indel_end = indel_start; //indel_start == first inserted/deleted base
				if (t[3].length == 1 && t[4].length > 1) is_deletion = false;
				if (is_deletion) indel_end = indel_start + t[3].length - 2; //indel_end == last deleted base
				else indel_end = indel_start + t[4].length - 2; //indel_end == last inserted base if aligned as soft clip
				if (flt_lv == 1) { //exact match btw indel and merged window
					if (indel_start == min_mg_start && indel_end == max_mg_end) flt = flt_dmg = true;
				}
				else if (flt_lv == 2) { //allow up to 2bp difference btw indel and merged window
					if (indel_start >= min_mg_start - 2 && indel_start <= min_mg_start + 2 && indel_end >= max_mg_end - 2 && indel_end <= max_mg_end + 2) flt = flt_dmg = true;
				}
				else if (flt_lv == 3) { //any overlap btw indel and merged window (default)
					if (indel_start <= max_mg_end && (indel_start >= min_mg_start || indel_end >= min_mg_start)) flt = flt_dmg = true;
				}
			} 
			if (cell[cell_id] == null) {
				cell[cell_id] = { flt:flt, dp:dp, ad:ad, adf:adf, adr:adr, lt:lt, flt_dmg:flt_dmg, count_ltpos:count_ltpos, uniq_ltpos:uniq_ltpos, min_mg_start:min_mg_start, max_mg_end:max_mg_end };
			} else {
				var c = cell[cell_id];
				if (flt) c.flt = flt;
				if (flt_dmg) c.flt_dmg = true;
				if (c.lt < lt) c.lt = lt;
				c.dp = 0;
				for (var j = 0; j < ad.length; ++j) {
					if (c.adf[j] > adf[j]) c.adf[j] = adf[j];
					if (c.adr[j] > adr[j]) c.adr[j] = adr[j];
					c.ad[j] = c.adf[j] + c.adr[j];
					c.dp += c.ad[j];
				}
				c.count_ltpos = count_ltpos;
				c.uniq_ltpos = uniq_ltpos;
				c.min_mg_start = min_mg_start;
				c.max_mg_end = max_mg_end;
			}
		}
	}

	// only consider beallelic sites for calling
	var alt = t[4].split(",");
	if (alt.length != 1 || (alt[0].length == 1 && t[3].length == 1)) //only keep biallelic indel sites
		flt_bulks = flt_indel = true;

	// test het in the bulk(s) //hom = ALT-hom in bulk
	var all_het = true, all_hom = true, all_good_alt = true;
	for (var i = 0; i < bulk.length; ++i) {
		var b = bulk[i];
		b.het = b.hom = false;
		if (b.adf[0] > 0 && b.adf[1] > 0 && b.adr[0] > 0 && b.adr[1] > 0 && b.ad[0] >= min_het_dp_bulk && b.ad[1] >= min_het_dp_bulk) {
			if (b.ad[0] >= b.dp * min_het_ab_bulk && b.ad[1] >= b.dp * min_het_ab_bulk)
				b.het = true;
		}
		if (!b.het && b.ad[1] >= min_het_dp_bulk && b.adf[1] > 0 && b.adr[1] > 0 && b.ad[0] <= max_alt_dp_bulk)
			b.hom = true;
		if (b.ad[1] < min_het_dp_bulk) all_good_alt = false;
		if (!b.het) all_het = false;
		if (!b.hom) all_hom = false;
		if (b.dp < min_dp_bulk)
			flt_bulks = true;
	}

	// output differences in bulk
	if (n_bulk > 1 && !is_hap_bulk) {
		var bulk_diff = false, n_bulk_ref = 0, n_bulk_alt = 0;
		for (var i = 0; i < bulk.length; ++i) {
			var b = bulk[i];
			if (b.ad[1] == 0) ++n_bulk_ref;
			else if (b.ad[1] >= min_dp_alt_cell && b.adf[1] >= min_dp_alt_strand_cell && b.adr[1] >= min_dp_alt_strand_cell)
				++n_bulk_alt;
		}
		if (n_bulk_ref > 0 && n_bulk_alt > 0) {
			var ad = [];
			for (var i = 0; i < bulk.length; ++i)
				ad.push(bulk[i].adf[1] + ':' + bulk[i].adr[1]);

			while (last_bulk.length && (last_bulk[0].ctg != t[0] || last_bulk[0].pos + flt_win < t[1])) {
				var x = last_bulk.shift();
				if (!x.flt) print('BV', x.data);
			}
			var flt_this = flt_bulks;
			if (var_map && var_map.get(t[0] + ':' + t[1]) != null)
				flt_this = true;
			for (var j = 0; j < last_bulk.length; ++j) {
				flt_this = true;
				last_bulk[j].flt = true;
			}
			last_bulk.push({ flt:flt_this, ctg:t[0], pos:t[1], data:[t[0], t[1], t[3], t[4], ad.join("\t")].join("\t") });
		}
	}

	// count ADO
	if (is_hap_bulk && all_hom && !flt_bulks) { //for hap bulk, only count alt allele dropped 
		++n_hom_bulk;
		for (var j = 0; j < cell.length; ++j)
			if (cell[j].flt || cell[j].ad[1] < min_joint_cell)
				++cell_meta[j].ado[1];
	}
	if (!is_hap_bulk && all_het && !flt_bulks) {
		++n_het_bulk;
		for (var j = 0; j < cell.length; ++j) {
			if (cell[j].flt || cell[j].ad[0] < min_joint_cell) ++cell_meta[j].ado[0]; // ref allele dropped
			if (cell[j].flt || cell[j].ad[1] < min_joint_cell) ++cell_meta[j].ado[1]; // alt allele dropped
		}
	}

	// test if ALT is callable and count FN
	var n_joint_alt = 0;
	var alt_detected = false;
	for (var i = 0; i < cell.length; ++i) {
		var c = cell[i];
		// If a cell is haploid and it has ref alleles, c.flt will be true. The conditions below work with haploid cells.
		c.alt = (!c.flt && c.ad[1] >= min_dp_alt_cell && c.adf[1] >= min_dp_alt_strand_cell && c.adr[1] >= min_dp_alt_strand_cell && c.ad[1] >= c.dp * min_ab_cell && c.ad[0] <= max_dp_ref_cell);
		c.joint_alt = (!c.flt && c.ad[1] >= min_joint_cell && c.adf[1] >= min_joint_strand_cell && c.adr[1] >= min_joint_strand_cell);
		if (c.joint_alt) ++n_joint_alt;
		// count FN when bulk is ALT but cell is not
		if (!flt_bulks && !c.alt && ((is_hap_bulk && all_hom) || (!is_hap_bulk && all_het)))
			++cell_meta[i].fn; 
		if (!flt_bulks && c.alt)
			alt_detected = true;
		// whether to call a damage
		c.dmg = (!c.flt_dmg && c.ad[1] >= min_dp_dmg_strand && c.ad[0] >= min_dp_dmg_strand && c.adf[1] * c.adr[1] == 0 && c.adf[0] * c.adf[1] == 0 && c.adr[0] * c.adr[1] == 0);
		if (all_het && !flt_bulks) {
			if (!(!c.flt_dmg && c.adf[0] >= min_dp_dmg_strand && c.adr[0] >= min_dp_dmg_strand)) ++cell_meta[i].dmg_fn[0]; //bulk het, but cell doesn't have enough REF support, count dmg FN for REF
			if (!(!c.flt_dmg && c.adf[1] >= min_dp_dmg_strand && c.adr[1] >= min_dp_dmg_strand)) ++cell_meta[i].dmg_fn[1]; //bulk het, but cell doesn't have enough ALT support, count dmg FN for ALT
			if (cell_meta[i].ploidy > 1 && c.dmg) ++cell_meta[i].dmg_fp; // no dmg_fp of this kind for a haploid cell //bulk het, but cell dmg (ALT only on one strand), count dmg FP
		}
	}

	// count detected germline in bulk
	if (!flt_bulks && alt_detected && ((is_hap_bulk && all_hom) || (!is_hap_bulk && all_het)))
		++n_het_bulk_detected;

	// skip the highly unlikely scenario: all bulks have good ALT alleles. The site is not used for window filtering.
	if (all_good_alt) continue; // bulk shouldn't be het (REF/ALT) at somatic variant sites, hence requiring at least one bulk with REF/REF below

	// requiring at least one bulk to have good RefHom
	var n_bulk_ref = 0;
	for (var i = 0; i < bulk.length; ++i)
		if (bulk[i].ad[1] <= max_alt_dp_bulk)
			++n_bulk_ref;
	if (n_bulk_ref == 0) flt_indel = true; // flag the infavorable scenario: no bulks with good RefHom; this site may be used for window filtering later

	// print sites with conflicting strand information
	if (!flt_indel && !flt_bulks) {
		var tmp = [];
		for (var i = 0; i < cell.length; ++i) {
			var c = cell[i];
			if (!c.dmg) continue;
			tmp.push(cell_meta[i].name + ':' + c.adf.join(",") + ':' + c.adr.join(",") + ':' + c.count_ltpos + ':' + c.uniq_ltpos + ':' + c.min_mg_start + ':' + c.max_mg_end); 
			++cell_meta[i].dmg;
		}
		if (tmp.length > 0 && var_map && var_map.get(t[0] + ':' + t[1]) != null) tmp.length = 0;
		if (tmp.length > 0) {
			// only consider bulk without strand info (same output format as NV)
			var bulk_ad = [0, 0], bulk_alt = [];
			for (var i = 0; i < bulk.length; ++i)
				bulk_ad[0] += bulk[i].ad[0], bulk_ad[1] += bulk[i].ad[1], bulk_alt.push(bulk[i].ad[1]);
			if (bulk_ad[1] != 0) bulk_ad[1] = bulk_alt.join(":");
			print('DV', t[0], t[1], t[3], t[4], bulk_ad.join("\t"), tmp.length, tmp.join(","));
		}
	}

	// test Indel
	var cell_alt_f = 0, cell_alt_r = 0;
	for (var i = 0; i < cell.length; ++i) {
		if (cell[i].flt) continue;
		cell_alt_f += cell[i].adf[1];
		cell_alt_r += cell[i].adr[1];
	}
	if (cell_alt_f < min_dp_alt_strand_cell || cell_alt_r < min_dp_alt_strand_cell || cell_alt_f + cell_alt_r < min_dp_alt_cell) // too few ALT reads in cell(s)
		flt_indel = true;

	// filter by window & print
	while (last.length && (last[0].ctg != t[0] || last[0].pos + flt_win < t[1])) {
		var x = last.shift();
		if (show_flt || !x.flt) aggregate_calls(x, cell_meta, is_hap_bulk);
	}

	var flt_this = flt_indel;
	if (flt_bulks) flt_this = true;
	if (var_map && var_map.get(t[0] + ':' + t[1]) != null)
		flt_this = true;
	for (var j = 0; j < last.length; ++j) { // go through rest of the candidates (failed window filter) and if same cell have two adjacent variants, mark both as filtered. Usure why keeping them in last[].
		for (var i = 0; i < cell.length; ++i)
			if (cell[i].ad[1] > 0 && last[j].cell[i].ad[1] > 0)
				flt_this = last[j].flt = true;
	}

	last.push({ flt:flt_this, n_joint_alt:n_joint_alt, ctg:t[0], pos:t[1], bulk:bulk, cell:cell, ref:t[3], alt:t[4] });
}
while (last_bulk.length) {
	var x = last_bulk.shift();
	if (!x.flt) print('BV', x.data);
}
while (last.length) {
	var x = last.shift();
	if (show_flt || !x.flt) aggregate_calls(x, cell_meta, is_hap_bulk);
}

/***************************
 * Output final statistics *
 ***************************/

var indel = [], fnr = [], corr_indel = [], dmg = [], corr_dmg = [], ado = [], fnr_dmg = [], fpr_dmg = [];
for (var i = 0; i < cell_meta.length; ++i) {
	var c = cell_meta[i];
	ado[i] = is_hap_bulk? c.ado[1] / n_hom_bulk : c.ploidy == 1? 2. * c.ado[1] / n_het_bulk - 1. : c.ado[1] / n_het_bulk;
	indel[i] = c.indel;
	fnr[i] = is_hap_bulk? c.fn / n_hom_bulk : c.ploidy == 1? 2. * c.fn / n_het_bulk - 1. : c.fn / n_het_bulk;
	corr_indel[i] = indel[i] / (1.0 - fnr[i]);
	dmg[i] = c.dmg;
	if (!is_hap_bulk) {
		fpr_dmg[i] = c.dmg_fp / n_het_bulk;
		fnr_dmg[i] = c.ploidy == 1? 2. * c.dmg_fn[1] / n_het_bulk - 1. : c.dmg_fn[1] / n_het_bulk;
		corr_dmg[i] = ((dmg[i] - corr_indel[i] * fpr_dmg[i]) / (1.0 - fnr[i])).toFixed(2);
		fnr_dmg[i] = fnr_dmg[i].toFixed(4);
		fpr_dmg[i] = fpr_dmg[i].toFixed(4);
	}
	fnr[i] = fnr[i].toFixed(4);
	corr_indel[i] = corr_indel[i].toFixed(2);
}
print('NN', indel.join("\t"));
print('NR', fnr.join("\t"));
print('NC', corr_indel.join("\t"));
print('DN', dmg.join("\t"));
if (!is_hap_bulk) {
	print('DP', fpr_dmg.join("\t"));
	print('DR', fnr_dmg.join("\t"));
	print('DC', corr_dmg.join("\t"));
}

// output "multi-alignment"
for (var i = 0; i < cell_meta.length; ++i)
	print('NA', cell_meta[i].name, ado[i].toFixed(4), cell_meta[i].calls.join(""));

// output sensitivity (binary or sum of TPR)
print('Sensitivity_binary', n_het_bulk_detected, n_het_bulk, (n_het_bulk_detected / n_het_bulk).toFixed(4));
var sum_fnr = fnr.reduce(function(a, b) { return parseFloat(a) + parseFloat(b); }, 0)
print('Sensitivity_FNR', fnr.length, sum_fnr, (fnr.length - sum_fnr).toFixed(4));

/********
 * Free *
 ********/

if (var_map != null) var_map.destroy();
buf.destroy();
file.close();
