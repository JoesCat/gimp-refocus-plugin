<?xml version='1.0'?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                version='1.0'>

<xsl:output method="xml"
  doctype-public="-//OASIS//DTD DocBook XML V4.1.2//EN"
  doctype-system="docbook/xml-dtd-4.1.2-1.0-8/docbookx.dtd"/>

<xsl:key name="id" match="*" use="@id"/>

<xsl:template match="/">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="*|@*">
  <xsl:copy>
  <xsl:apply-templates select="@*|node()"/>
</xsl:copy>
</xsl:template>

<xsl:template match = "xref">
  <xsl:value-of select="id(./@linkend)/@xreflabel"/>
</xsl:template>

<xsl:template match = "ulink">
  <xsl:value-of select="./@url"/>
</xsl:template>
</xsl:stylesheet>
