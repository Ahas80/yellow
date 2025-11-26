BEGIN{FS=OFS="\t"}
/^#/ {print; next}
{
  # drop ultra-short CDS (<270 nt == 90 aa)
  if ($3=="CDS") { len=$5-$4+1; if (len<270) next }
  attr=$9
  id=""
  if (match(attr,/ID=[^;]+/))         id=substr(attr,RSTART+3,RLENGTH-3)
  if (id=="") {
    if (match(attr,/locus_tag=[^;]+/)) {
      lt=substr(attr,RSTART+10,RLENGTH-10)
      if (lt!="") attr="ID="lt";"attr
      id=lt
    }
  }
  # if still no ID, keep line as-is
  if (id=="") { $9=attr; print; next }
  # de-duplicate IDs within file
  count[id]++
  if (count[id]>1) {
    new_id=id"__dup"count[id]
    gsub("ID="id, "ID="new_id, attr)
  }
  $9=attr
  print
}
