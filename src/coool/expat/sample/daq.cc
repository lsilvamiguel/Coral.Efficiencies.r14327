/* This is simple demonstration of how to use expat. This program
reads an XML document from standard input and writes a line with the
name of each element to standard output indenting child elements by
one tab stop more than their parent element. */

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "xmlparse.h"

char *prefix(int depth)
{
  static char s[1111];
  s[0]=0;
  while( depth-- > 0 )
    strcat(s,"        ");
  return s;
}

/* XML_StartElementHandler */
void startElement(void *userData, const char *name, const char **atts)
{
  int i;
  int *depthPtr = (int*)userData;

  printf("%s%s\n",prefix(*(int*)userData),name);
  for( int i=0; atts[i]; i++ )
    printf("  %s%s\n",prefix(*(int*)userData),atts[i]);

  *depthPtr += 1;
}

/* XML_EndElementHandler */
void endElement(void *userData, const char *name)
{
  int *depthPtr = (int*)userData;
  *depthPtr -= 1;
}

/* XML_CharacterDataHandler */
void characters_handler(void *userData,const XML_Char *s,int len)
{
  printf("%s\"%s\"\n",prefix(*(int*)userData),s);
}

/* XML_ProcessingInstructionHandler */
void instruction_handler(void *userData,const XML_Char *target,const XML_Char *data)
{
  printf("%sInstruction:  target=%s   data=%s\n",prefix(*(int*)userData),target,data);
}

/* XML_DefaultHandler */
void default_handler(void *userData,const XML_Char *s,int len)
{
  char *p,ss[len+1];
  strncpy(ss,s,len);
  ss[len]=0;
  for( p=ss; p<ss+len; p++ )
  {
    if( *p==' ' || *p=='	' )  // space or tab
      continue;
    else
    if( *p=='\n' )
    {
     /*
      *       if( (p-ss+1)!=len )
      *       {
      *         printf("(p-ss+1)=%d len=%d  \\n=%d  c=%d\n",p-ss+1,len,int('\n'),p[1]);
      *       }
      */
      return;
    }
    else
      break;
  }

  if( p==ss+len )
    return; // all spaces

  if( strchr(ss,'\n')==ss+len )
    return;
  printf("%s------------------------\n",prefix(*(int*)userData));
  printf("%s[%3d]\"%s\"\n",prefix(*(int*)userData),len,ss);
  printf("%s------------------------\n",prefix(*(int*)userData));
}

/* XML_CommentHandler */
void comment_handler(void *userData, const XML_Char *data)
{
}





int main()
{
  char buf[BUFSIZ];
  XML_Parser parser = XML_ParserCreate(NULL);
  int done;
  int depth = 0;
  XML_SetUserData                       (parser, &depth);
  XML_SetElementHandler                 (parser, startElement, endElement);
//XML_SetCharacterDataHandler           (parser, characters_handler);
  XML_SetProcessingInstructionHandler   (parser, instruction_handler);
  XML_SetDefaultHandler                 (parser, default_handler);
  XML_SetCommentHandler                 (parser, comment_handler);
  do {
    size_t len = fread(buf, 1, sizeof(buf), stdin);
    done = len < sizeof(buf);
    if (!XML_Parse(parser, buf, len, done)) {
      fprintf(stderr,
	      "%s at line %d\n",
	      XML_ErrorString(XML_GetErrorCode(parser)),
	      XML_GetCurrentLineNumber(parser));
      return 1;
    }
  } while (!done);
  XML_ParserFree(parser);
  return 0;
}
