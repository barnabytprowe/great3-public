#coding: utf-8
from django import template
from leaderboard.models import PLACEHOLDER_ERROR

register = template.Library()

@register.filter
def with_1e3_error(value, arg):
	if float(value)==PLACEHOLDER_ERROR:
		return ""
	return "%.3g ± %.3g" % (1e3*value, 1e3*arg)

@register.filter
def with_1e4_error(value, arg):
	if float(value)==PLACEHOLDER_ERROR:
		return ""
	return "%.3g ± %.3g" % (1e4*value, 1e4*arg)


@register.filter
def split_lines(value):
	value = str(value)
	if len(value)>16:
		value = value[:16] + "\n" + value[16:]
		print value
	return value


@register.filter
def shorten_name(value):
	return str(value).replace("_", " ")
