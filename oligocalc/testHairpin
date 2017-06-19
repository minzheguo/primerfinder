function MakeComplement(theSequence, isDNA) {
	var returnString="";
	var i;
	var temp;
	for( i=theSequence.length-1; i>=0; i--) {
		temp=theSequence.charAt(i);
		switch (temp) {
			case "A" :
				if (isDNA) {
					temp="T";
				} else {
					temp="U";
				}
				break;
			case "T" :
				temp="A";
				break;
			case "U" :
				temp="A";
				break;
			case "G" :
				temp="C";
				break;
			case "C" :
				temp="G";
				break;
			case "M" :
				temp="K";
				break;
			case "K" :
				temp="M";
				break;
			case "R" :
				temp="Y";
				break;
			case "Y" :
				temp="R";
				break;
			case "W" :
				temp="W";
				break;
			case "S" :
				temp="S";
				break;
			case "V" :
				temp="B";
				break;
			case "B" :
				temp="V";
				break;
			case "H" :
				temp="D";
				break;
			case "D" :
				temp="H";
				break;
			default : break;
		}
		returnString=returnString+temp;
	}
	return returnString;
}

function getIndexOf(seq, subSeq, startIndex, minMatch, broadMatch=false)
{
// look for subSeq in seq
/* returns an array where
	theResult[0] is the index of the first match of subseq that is of at least length minMatch in seq
	theResult[1] is the length of the match
*/
	var theResult= new Array(2);
	theResult[0]=-1;
	theResult[1]=-1;
	if (!broadMatch) {
		for(var k=minMatch; k<=subSeq.length; k++) {
			// can replace this with seq.search for GREP capabilities
			theMatch=seq.indexOf(subSeq.substring(0,k),startIndex);
			if (theMatch < 0) {
				break;
			} else {
				theResult[0]= theMatch;
				theResult[1] = k;
				
			}
		}
	} else {
		for(var i=startIndex; i<seq.length; i++) {
			if(isBaseEqual(seq.charAt(i),subSeq.charAt(0))) {
		 		for(j=0; j<subSeq.length; j++) {
					if(!isBaseEqual(seq.charAt(i+j),subSeq.charAt(j))) {
						break;
					} else if (j >= minMatch-1) {
						theResult[0]= theMatch;
						theResult[1] = k;
					}
				}	
				if (j==subSeq.length) {
						theResult[0]= theMatch;
						theResult[1] = k;
				}
			}
		}
	}
	
	return theResult;
}

function DoHairpinArrayInsert(a,b,c,d,results)
{
	var arrayCount=results.length;
	if (a >= c || a >= b || c >= d || b>=c) {
		return results;
	}
	for (var i=0;i<arrayCount;i++) {
		if (results[i][0]<=a && results[i][1]>=b && results[i][2]<=c && results[i][3]>=d)
			return results;
		if (results[i][0]>=a && results[i][1]<=b && results[i][2]>=c && results[i][3]<=d) {
			results[i][0]=a;
			results[i][1]=b;
			results[i][2]=c;
			results[i][3]=d;
			return results;
		}
	}
	results[arrayCount]=new Array(4);
	results[arrayCount][0]=a;
	results[arrayCount][1]=b;
	results[arrayCount][2]=c;
	results[arrayCount][3]=d;
	return results;
}

var theFullSequence='AGAGCCAAAGAAACTGGCAA';
var minHairpinLength=4;
var bubbleSize=3;

var theFullComplement=MakeComplement(theFullSequence, 1);

var theResults = new Array();


	var theResult;
	var count;
	var compPos;
	var seqPos;
	var maxSeqLength=Math.abs(theFullSequence.length/2)-bubbleSize; 
	
	var maxMatch=0;
	for (compPos=0; compPos<theFullComplement.length-2*minHairpinLength; compPos++) {
		maxMatch=0;
		for (seqPos=0; seqPos<theFullSequence.length-maxSeqLength; seqPos++) {
					             theResult=getIndexOf(theFullSequence.substring(0,seqPos+maxSeqLength),theFullComplement.substring(compPos,theFullComplement.length), seqPos, minHairpinLength);
			if (theResult[0] > -1) {
				 alert(theResult)
				 theResults=DoHairpinArrayInsert(theResult[0],theResult[0]+theResult[1]-1,theFullSequence.length-compPos-theResult[1],theFullSequence.length-compPos-1,theResults);

				 if (theResult[1] > maxMatch) maxMatch=theResult[1];
				seqPos=theResult[0]+theResult[1]-minHairpinLength;  // move forward to guarantee nothing else is found that is a reasonable match
				if (seqPos+minHairpinLength>=maxSeqLength) {
					compPos+=maxMatch-minHairpinLength; // move compPos forward to stop identical checks if long match was found!
					break; // we have moved far enough on the primer to guarentee we have everything  -further would give us the reverse match
				}
			} else {
				if (maxMatch > minHairpinLength) compPos+=maxMatch-minHairpinLength; 
				break;  //not found in the rest of the sequence!
			}
		}
	}
	
alert(theResults)
