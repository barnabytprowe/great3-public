
function adjustFieldVisibility(type, listener){
	var methodField = document.getElementById("id_method");
	var chosenMethod = methodField.options[methodField.selectedIndex].value;
	if (chosenMethod=='None'){
		document.getElementById('id_new_method_name').parentNode.style.visibility='visible';
 	} else {
		document.getElementById('id_new_method_name').parentNode.style.visibility='hidden';
	}
}

var methodFieldThing = document.getElementById("id_method");
if (methodFieldThing.addEventListener) {
	methodFieldThing.addEventListener("change", adjustFieldVisibility, false);
} else {
	methodFieldThing.attachEvent('onchange', adjustFieldVisibility);
}  

window.onload = function(e){ 
	adjustFieldVisibility(0,0);
}
